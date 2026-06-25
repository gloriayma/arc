"""CoreWithRNA: subclass of NTv3's Core that accepts an optional RNA bias
at the post-conv / pre-transformer bottleneck (B, L/128, embed_dim).

The forward body is copy-pasted verbatim from
    ntv3/ntv3_base_model/modeling_ntv3_pretrained.py :: Core.forward
with one additional block marked `# >>> sc_nona >>>` adding `rna_bias` at the
L/128 bottleneck. Keep it in sync with upstream if NTv3 source changes.
"""
from __future__ import annotations

import torch
from torch.nn import functional as F

# NTv3 source is on sys.path via model/__init__.py.
from configuration_ntv3_pretrained import Ntv3PreTrainedConfig
from modeling_ntv3_pretrained import (
    Core,
    _autocast_to,
    _dtype_from_str,
)


class CoreWithRNA(Core):
    """NTv3 Core with an optional `rna_bias` added at the L/128 bottleneck.

    `rna_bias` shape must be `(B, L/128, config.embed_dim)`. When `rna_bias`
    is None or zeros, this module is bit-exact-equivalent to stock NTv3 Core.
    """

    def forward(
        self,
        input_ids: torch.LongTensor | None = None,
        inputs_embeds: torch.FloatTensor | None = None,
        rna_bias: torch.Tensor | None = None,
        output_hidden_states: bool = False,
        output_attentions: bool = False,
    ) -> dict[str, torch.Tensor | list[torch.Tensor]]:
        assert (input_ids is None) != (inputs_embeds is None), \
            "You must specify exactly one of input_ids or inputs_embeds"

        device_type = (
            input_ids.device.type
            if input_ids is not None
            else inputs_embeds.device.type  # type: ignore
        )
        hidden_states: list[torch.Tensor] = []
        attentions: list[torch.Tensor] = []
        outs: dict[str, torch.Tensor | list[torch.Tensor]] = {}

        # Embedding
        if inputs_embeds is None:
            x = self.embed_layer(input_ids)
        else:
            x = inputs_embeds
        emb_compute = _dtype_from_str(
            getattr(self.config, "embedding_compute_dtype", "float32")
        )
        x = x.to(emb_compute)

        # Stem
        with _autocast_to(
            device_type, getattr(self.config, "stem_compute_dtype", "bfloat16")
        ):
            x = self.stem(x.permute(0, 2, 1))

        # Down conv tower
        with _autocast_to(
            device_type,
            getattr(self.config, "down_convolution_compute_dtype", "bfloat16"),
        ):
            residuals: list[torch.Tensor] = []
            for block in self.conv_tower_blocks:
                residuals.append(x)
                y = block.conv(x)
                y = block.res_conv(y)
                if output_hidden_states:
                    hidden_states.append(y.permute(0, 2, 1))
                x = block.avg_pool(y)

        # Transformer tower
        x = x.permute(0, 2, 1)  # (B, L/128, C)
        # >>> sc_nona >>> RNA bias injection at the post-conv / pre-transformer bottleneck
        if rna_bias is not None:
            x = x + rna_bias.to(x.dtype)
        # <<< sc_nona <<<
        with _autocast_to(
            device_type,
            getattr(self.config, "transformer_qkvo_compute_dtype", "bfloat16"),
        ):
            for i, layer in enumerate(self.transformer_blocks):
                out = layer(x, attention_mask=None, attention_weight_bias=None)
                x = out["embeddings"]
                if output_hidden_states:
                    hidden_states.append(x)
                if output_attentions:
                    attentions.append(out["attention_weights"])
                if (i + 1) in self.config.embeddings_layers_to_save:
                    outs[f"embeddings_{i + 1}"] = x
                if (i + 1) in self._attention_layers_to_save:
                    for m in self._attention_maps_per_layer_to_save[i + 1]:
                        outs[f"attention_map_layer_{i + 1}_number_{m}"] = out[
                            "attention_weights"
                        ][:, m + 1]

        # Deconv tower
        x = x.permute(0, 2, 1)  # (B, C, L/128)
        with _autocast_to(
            device_type,
            getattr(self.config, "up_convolution_compute_dtype", "bfloat16"),
        ):
            for i, block in enumerate(self.deconv_tower_blocks):
                r = residuals[-(i + 1)]
                y = block(x)
                if self.config.use_skip_connection:
                    y = y + r
                x = y
                if output_hidden_states:
                    hidden_states.append(x.permute(0, 2, 1))
                if (i + 1) in self.config.deconv_layers_to_save:
                    outs[f"embeddings_deconv_{i + 1}"] = x

        # Head
        y = x.permute(0, 2, 1)  # (B, L, C)
        with _autocast_to(
            device_type, getattr(self.config, "lmhead_compute_dtype", "float32")
        ):
            y = F.gelu(y, approximate="tanh")
            for hl in self.lm_head["hidden_layers"]:
                y = F.gelu(hl(y), approximate="tanh")
            y = y.to(self.lm_head["head"].weight.dtype)
            logits = self.lm_head["head"](y)

        outs["logits"] = logits
        if output_hidden_states:
            outs["hidden_states"] = hidden_states
        if output_attentions:
            outs["attentions"] = attentions
        return outs


def load_pretrained(model_name: str, *, token: str | None = None) -> CoreWithRNA:
    """Instantiate CoreWithRNA and load weights from an NTv3 HF checkpoint.

    Strips the `core.` prefix on state_dict keys because the saved NTv3 model
    is `NTv3PreTrained.core`, not `Core` directly.
    """
    from huggingface_hub import hf_hub_download
    from safetensors.torch import load_file

    cfg = Ntv3PreTrainedConfig.from_pretrained(model_name, token=token)
    model = CoreWithRNA(cfg)
    weights_path = hf_hub_download(model_name, "model.safetensors", token=token)
    sd_full = load_file(weights_path)
    sd = {k.removeprefix("core."): v for k, v in sd_full.items() if k.startswith("core.")}
    # strict=False lets the lazily-built rotary cache buffers be (re)created on
    # first forward without producing load_state_dict errors.
    model.load_state_dict(sd, strict=False)
    return model
