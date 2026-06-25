"""RNAEncoder built out of NTv3 components wherever possible.

NTv3 parts reused (imported from ntv3_base_model/modeling_ntv3_pretrained.py):
    - LayerNormFP32                (always fp32 internally, torch.compile-friendly)
    - SelfAttentionBlock           (pre-norm + MHA + GLU-SiLU FFN, with optional RoPE)
    - MultiHeadAttention + RotaryEmbedding   (used internally by SelfAttentionBlock)
    - F.gelu(..., approximate="tanh")        (matches NTv3's stem/lm_head activation)

So the only "new" mechanisms specific to this encoder are:
    - mean-pool over L (NTv3 uses strided convs to downsample; we match the L/128 ratio)
    - log1p of the binned values
    - per-(cell, bin) value MLP (scalar -> d_rna)
    - learned mask embedding for masked entries
    - axial structure: run a SelfAttentionBlock over L (with RoPE), then over D (no RoPE)
    - learned-query collapse of the D axis
    - zero-init projection to d_ntv -> rna_bias
"""
from __future__ import annotations

import math

import torch
from torch import nn
from torch.nn import functional as F

# NTv3 source is on sys.path via model/__init__.py.
from modeling_ntv3_pretrained import (
    LayerNormFP32,
    RotaryEmbeddingConfig,
    SelfAttentionBlock,
)


class ValueMLP(nn.Module):
    """Shared (scalar value) -> (d_rna vector) map. NTv3-style GELU(tanh).
    Basically just an embedding layer for the per-cell per-"base" scalar coverage signal."""


    def __init__(self, d_rna: int):
        super().__init__()
        self.fc1 = nn.Linear(1, d_rna)
        self.fc2 = nn.Linear(d_rna, d_rna)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = F.gelu(self.fc1(x), approximate="tanh")
        return self.fc2(x)


class AxialAttentionBlock(nn.Module):
    """One axial block: SelfAttentionBlock over L axis (with RoPE), then over D (no RoPE).

    SelfAttentionBlock includes its own residuals + pre-LayerNormFP32 + GLU-SiLU FFN,
    so this wrapper just handles the reshape/permute around each call.
    """

    def __init__(
        self,
        d_rna: int,
        num_heads: int = 4,
        ffn_embed_dim: int | None = None,
        layer_norm_eps: float = 1e-5,
    ):
        super().__init__()
        ffn_embed_dim = ffn_embed_dim or 4 * d_rna
        # L-axis attention uses RoPE: positions are ordered, distances meaningful.
        self.attn_l = SelfAttentionBlock(
            num_heads=num_heads,
            embed_dim=d_rna,
            ffn_embed_dim=ffn_embed_dim,
            add_bias_kv=False,
            add_bias_fnn=False,
            ffn_activation_name="swish",
            use_glu_in_ffn=True,
            layer_norm_eps=layer_norm_eps,
            pre_layer_norm=True,
            rotary_embedding_config=RotaryEmbeddingConfig(rescaling_factor=None),
        )
        # D-axis attention: cells are exchangeable -> no positional encoding.
        self.attn_d = SelfAttentionBlock(
            num_heads=num_heads,
            embed_dim=d_rna,
            ffn_embed_dim=ffn_embed_dim,
            add_bias_kv=False,
            add_bias_fnn=False,
            ffn_activation_name="swish",
            use_glu_in_ffn=True,
            layer_norm_eps=layer_norm_eps,
            pre_layer_norm=True,
            rotary_embedding_config=None,
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # x: (B, L, D, d_rna)
        B, L, D, d_rna = x.shape

        # L-axis: each (b, d) is a sequence of length L.
        h = x.permute(0, 2, 1, 3).reshape(B * D, L, d_rna)
        h = self.attn_l(h)["embeddings"]
        x = h.reshape(B, D, L, d_rna).permute(0, 2, 1, 3)

        # D-axis: each (b, l) is a sequence of length D.
        h = x.reshape(B * L, D, d_rna)
        h = self.attn_d(h)["embeddings"]
        x = h.reshape(B, L, D, d_rna)
        return x


class RNAEncoder(nn.Module):
    def __init__(
        self,
        d_rna: int,
        d_ntv: int,
        num_axial_layers: int = 4,
        num_heads: int = 4,
        ffn_embed_dim: int | None = None,
        bin_size: int = 128,
        layer_norm_eps: float = 1e-5,
    ):
        super().__init__()
        self.d_rna = d_rna
        self.d_ntv = d_ntv
        self.bin_size = bin_size

        self.value_mlp = ValueMLP(d_rna)

        # Learned mask embedding for masked (cell, bin) entries.
        self.mask_emb = nn.Parameter(torch.empty(d_rna))
        nn.init.normal_(self.mask_emb, std=0.02)

        self.blocks = nn.ModuleList(
            [
                AxialAttentionBlock(
                    d_rna=d_rna,
                    num_heads=num_heads,
                    ffn_embed_dim=ffn_embed_dim,
                    layer_norm_eps=layer_norm_eps,
                )
                for _ in range(num_axial_layers)
            ]
        )
        self.norm_out = LayerNormFP32(d_rna, eps=layer_norm_eps)

        # Learned query for the D-axis collapse.
        self.collapse_query = nn.Parameter(torch.empty(d_rna))
        nn.init.normal_(self.collapse_query, std=0.02)

        # Zero-init projection to NTv3's d_ntv so day-0 rna_bias == 0.
        self.to_ntv = nn.Linear(d_rna, d_ntv)
        nn.init.zeros_(self.to_ntv.weight)
        nn.init.zeros_(self.to_ntv.bias)

    def forward(
        self,
        rna: torch.Tensor,
        rna_mask: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            rna:      (B, L, D) raw counts / coverage. L must be divisible by bin_size.
            rna_mask: (B, L/bin, D) bool, True at masked (cell, bin) entries to be predicted.
                      None = no masking applied.

        Returns:
            rna_hidden: (B, L/bin, D, d_rna)
            rna_bias:   (B, L/bin, d_ntv)
        """
        B, L, D = rna.shape
        assert L % self.bin_size == 0, f"L={L} not divisible by bin_size={self.bin_size}"
        L_bin = L // self.bin_size

        # Mean-pool L -> L_bin, log1p (counts are heavy-tailed).
        binned = rna.transpose(1, 2)  # (B, D, L)
        binned = F.avg_pool1d(binned, kernel_size=self.bin_size)  # (B, D, L_bin)
        binned = binned.transpose(1, 2)  # (B, L_bin, D)
        binned = torch.log1p(binned.clamp(min=0))

        # Per-(cell, bin) value MLP.
        x = self.value_mlp(binned.unsqueeze(-1))  # (B, L_bin, D, d_rna)

        # Masked entries get the mask embedding.
        if rna_mask is not None:
            assert rna_mask.shape == (B, L_bin, D), (
                f"rna_mask shape {tuple(rna_mask.shape)} != ({B}, {L_bin}, {D})"
            )
            mask_b = rna_mask.unsqueeze(-1)  # (B, L_bin, D, 1)
            x = torch.where(mask_b, self.mask_emb.to(x.dtype).expand_as(x), x)

        # No absolute positional embedding: RoPE inside the L-axis attention handles it.

        for blk in self.blocks:
            x = blk(x)
        x = self.norm_out(x)
        rna_hidden = x  # (B, L_bin, D, d_rna)

        # Collapse D via dot-product attention with a single learned query.
        scores = (x * self.collapse_query).sum(dim=-1) / math.sqrt(self.d_rna)  # (B, L_bin, D)
        weights = F.softmax(scores, dim=-1).unsqueeze(-1)  # (B, L_bin, D, 1)
        collapsed = (weights * x).sum(dim=-2)  # (B, L_bin, d_rna)

        rna_bias = self.to_ntv(collapsed)  # (B, L_bin, d_ntv)
        return rna_hidden, rna_bias
