"""Permutation-equivariant RNA-MLM head.

Reconstructs the per-(cell, bin) RNA value from:
    rna_hidden: (B, L_bin, D, d_rna)   uncollapsed encoder output (carries cell identity at each bin)
    dna_repr:   (B, L_bin, d_ntv)       DNA representation at the same bin

The same MLP is applied to every (b, l, d) tuple, so the output is identical
under any permutation of the D axis. This is what makes the head correct for
exchangeable-cell training.

Predictions are in **log1p space** (matching the RNAEncoder's input
preprocessing), so the training loss is plain MSE against log1p(target).
NTv3-style components are reused: LayerNormFP32, F.gelu(approximate="tanh").

TODO: prediction is at bin-level (L/128). Consider upgrading to 1bp resolution
later -- mirroring NTv3's bigwig head, which predicts after the deconv tower at
single-nucleotide resolution. Would need either (a) a parallel deconv tower for
the RNA stream, or (b) repeat-interleave the RNA hidden to L and pair with the
deconv'd DNA. Bin-level is fine for v1 (gene/peak granularity); 1bp would matter
for splice-site / single-base-resolution work.
"""
from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F

# NTv3 source is on sys.path via model/__init__.py.
from modeling_ntv3_pretrained import LayerNormFP32


class RNAMlmHead(nn.Module):
    def __init__(
        self,
        d_rna: int,
        d_ntv: int,
        hidden_dim: int | None = None,
        layer_norm_eps: float = 1e-5,
    ):
        super().__init__()
        hidden_dim = hidden_dim or d_rna

        # Project DNA bottleneck to the same dim as cell rows for concat.
        self.dna_proj = nn.Linear(d_ntv, d_rna)
        # Norm + MLP.
        self.norm = LayerNormFP32(2 * d_rna, eps=layer_norm_eps)
        self.fc1 = nn.Linear(2 * d_rna, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, 1)

    def forward(self, rna_hidden: torch.Tensor, dna_repr: torch.Tensor) -> torch.Tensor:
        """
        Args:
            rna_hidden: (B, L_bin, D, d_rna)
            dna_repr:   (B, L_bin, d_ntv)
        Returns:
            (B, L_bin, D) scalar predictions in log1p space.
        """
        B, L_bin, D, d_rna = rna_hidden.shape
        # Project DNA bottleneck and broadcast over D.
        dna = self.dna_proj(dna_repr)                       # (B, L_bin, d_rna)
        dna = dna.unsqueeze(2).expand(B, L_bin, D, d_rna)   # (B, L_bin, D, d_rna)

        x = torch.cat([rna_hidden, dna], dim=-1)            # (B, L_bin, D, 2*d_rna)
        x = self.norm(x)
        x = F.gelu(self.fc1(x), approximate="tanh")
        x = self.fc2(x).squeeze(-1)                         # (B, L_bin, D)
        return x
