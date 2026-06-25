"""NTv3WithRNA: composes CoreWithRNA + RNAEncoder + RNAMlmHead + DNARNAHead.

Forward: (dna_tokens, rna_matrix, rna_mask) -> {
    "dna_logits_stock":    (B, L, alphabet)     -- NTv3's lm_head applied to deconv'd DNA
    "dna_logits_with_rna": (B, L, alphabet)     -- DNARNAHead: DNA + upsampled rna_bias
    "rna_pred":            (B, L/bin, D)        -- RNAMlmHead per-cell continuous reconstruction
    "rna_bias":            (B, L/bin, d_ntv)    -- exposed for debugging / inspection
}

Day-0 (no training): rna_bias is exactly zero (RNAEncoder.to_ntv is zero-init), AND
DNARNAHead is initialized via from_pretrained_lm_head so its RNA-input weight slice
is zero. Result: both dna_logits_* match stock NTv3 bit-exactly.
"""
from __future__ import annotations

import torch
from torch import nn

# NTv3 source is on sys.path via model/__init__.py.
from configuration_ntv3_pretrained import Ntv3PreTrainedConfig

from model.core_with_rna import CoreWithRNA, load_pretrained
from model.dna_rna_head import DNARNAHead
from model.rna_encoder import RNAEncoder
from model.rna_mlm_head import RNAMlmHead


class NTv3WithRNA(nn.Module):
    def __init__(
        self,
        core: CoreWithRNA,
        d_rna: int = 64,
        num_axial_layers: int = 4,
        num_heads: int = 4,
    ):
        super().__init__()
        self.core = core
        cfg: Ntv3PreTrainedConfig = core.config
        self.config = cfg
        self.bin_size = 2 ** cfg.num_downsamples

        self.rna_encoder = RNAEncoder(
            d_rna=d_rna,
            d_ntv=cfg.embed_dim,
            num_axial_layers=num_axial_layers,
            num_heads=num_heads,
            bin_size=self.bin_size,
        )

        # Per-cell reconstruction head reads encoder rna_hidden + DNA at L/bin.
        self.rna_mlm_head = RNAMlmHead(d_rna=d_rna, d_ntv=cfg.embed_dim)

        # DNARNAHead initialized from stock lm_head, RNA slice zero -> day-0 identity.
        # d_rna_inject = embed_dim because we reuse rna_bias (already at d_ntv).
        self.dna_rna_head = DNARNAHead.from_pretrained_lm_head(
            config=cfg,
            d_rna_inject=cfg.embed_dim,
            pretrained_lm_head=core.lm_head,
        )

    @classmethod
    def from_pretrained(
        cls,
        model_name: str,
        *,
        d_rna: int = 64,
        token: str | None = None,
        num_axial_layers: int = 4,
        num_heads: int = 4,
    ) -> "NTv3WithRNA":
        core = load_pretrained(model_name, token=token)
        return cls(
            core=core,
            d_rna=d_rna,
            num_axial_layers=num_axial_layers,
            num_heads=num_heads,
        )

    def forward(
        self,
        dna_tokens: torch.Tensor,
        rna_matrix: torch.Tensor,
        rna_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        """
        Args:
            dna_tokens: (B, L) int token ids
            rna_matrix: (B, L, D) continuous per-base RNA coverage
            rna_mask:   (B, L/bin, D) bool or None, True at masked (cell, bin) entries
        """
        L = dna_tokens.shape[1]
        L_bin = L // self.bin_size

        rna_hidden, rna_bias = self.rna_encoder(rna_matrix, rna_mask=rna_mask)

        core_out = self.core(
            dna_tokens,
            rna_bias=rna_bias,
            output_hidden_states=True,
        )
        dna_logits_stock = core_out["logits"]                # (B, L, alphabet)
        hidden_states = core_out["hidden_states"]

        # DNA at L/bin = the last hidden state whose length axis equals L_bin
        # (= post-final-transformer-layer, pre-deconv). Robust to num_layers / num_downsamples.
        dna_at_bin = next(h for h in reversed(hidden_states) if h.shape[1] == L_bin)
        # DNA at full L = the last deconv output. Same tensor that stock lm_head consumed.
        dna_at_full = hidden_states[-1]

        dna_logits_with_rna = self.dna_rna_head(dna_at_full, rna_bias)
        rna_pred = self.rna_mlm_head(rna_hidden, dna_at_bin)

        return {
            "dna_logits_stock": dna_logits_stock,
            "dna_logits_with_rna": dna_logits_with_rna,
            "rna_pred": rna_pred,
            "rna_bias": rna_bias,
        }
