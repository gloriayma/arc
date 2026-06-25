"""DNARNAHead: NTv3-style DNA-MLM head that also reads the RNA stream.

Architecture mirrors NTv3's Core.lm_head exactly:
    input -> GELU(tanh) -> [num_hidden_layers_head x (Linear + GELU)] -> final Linear -> logits

The only difference is that the input at the first stage is the concatenation of
    dna_repr   : (B, L, embed_dim)         post-deconv DNA representation (1 bp resolution)
    rna_stream : (B, L/bin, d_rna_inject)  upsampled to (B, L, d_rna_inject) by nearest-repeat

Output: (B, L, alphabet_size) DNA-token logits, same as the stock head.

Use `DNARNAHead.from_pretrained_lm_head(config, d_rna_inject, stock_lm_head)` to
initialize from an NTv3 lm_head: the DNA-input weight slice is copied, the
RNA-input slice is zero-initialized, so day-0 logits == applying the stock head
to dna_repr alone (regardless of rna_stream).
"""
from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F


class DNARNAHead(nn.Module):
    def __init__(self, config, d_rna_inject: int):
        super().__init__()
        embed_dim = config.embed_dim
        final_in_dim = config.conv_init_embed_dim  # equals embed_dim in all NTv3 configs we've seen
        num_hidden = config.num_hidden_layers_head
        alphabet_size = config.alphabet_size

        self.bin_size = 2 ** config.num_downsamples
        self.d_rna_inject = d_rna_inject
        self.embed_dim = embed_dim

        if num_hidden > 0:
            hidden: list[nn.Linear] = []
            in_dim = embed_dim + d_rna_inject
            for _ in range(num_hidden):
                hidden.append(nn.Linear(in_dim, embed_dim))
                in_dim = embed_dim
            self.hidden_layers = nn.ModuleList(hidden)
            head_in_dim = final_in_dim
        else:
            # No hidden layers: the final Linear reads the concat directly.
            self.hidden_layers = nn.ModuleList()
            head_in_dim = final_in_dim + d_rna_inject

        self.head = nn.Linear(head_in_dim, alphabet_size)

    def forward(self, dna_repr: torch.Tensor, rna_stream: torch.Tensor) -> torch.Tensor:
        """
        Args:
            dna_repr:   (B, L, embed_dim) post-deconv DNA representation.
            rna_stream: (B, L/bin, d_rna_inject) collapsed RNA representation at L/bin resolution.
        Returns:
            (B, L, alphabet_size) DNA-token logits.
        """
        # Upsample RNA stream from L/bin to L by nearest-neighbor repeat.
        rna_up = rna_stream.repeat_interleave(self.bin_size, dim=1)
        rna_up = rna_up[:, : dna_repr.shape[1]]  # safety crop if shape mismatch

        x = torch.cat([dna_repr, rna_up], dim=-1)
        x = F.gelu(x, approximate="tanh")
        for hl in self.hidden_layers:
            x = F.gelu(hl(x), approximate="tanh")
        x = x.to(self.head.weight.dtype)
        return self.head(x)

    @classmethod
    def from_pretrained_lm_head(
        cls,
        config,
        d_rna_inject: int,
        pretrained_lm_head: nn.ModuleDict,
    ) -> "DNARNAHead":
        """Copy weights from a stock NTv3 lm_head into a new DNARNAHead.

        The DNA-input slice of each Linear's weight is copied verbatim; the
        RNA-input slice is zero-initialized. Result: day-0 output is identical
        to applying the stock head to dna_repr alone, regardless of rna_stream.
        """
        h = cls(config, d_rna_inject)
        embed_dim = config.embed_dim
        final_in_dim = config.conv_init_embed_dim

        stock_hidden = pretrained_lm_head["hidden_layers"]
        stock_head = pretrained_lm_head["head"]

        with torch.no_grad():
            for i, hl in enumerate(h.hidden_layers):
                stock = stock_hidden[i]
                hl.weight.zero_()
                if i == 0:
                    # First hidden layer's input is (embed_dim + d_rna_inject).
                    hl.weight[:, :embed_dim].copy_(stock.weight)
                else:
                    hl.weight.copy_(stock.weight)
                hl.bias.copy_(stock.bias)

            h.head.weight.zero_()
            if len(h.hidden_layers) == 0:
                # Final layer reads (final_in_dim + d_rna_inject).
                h.head.weight[:, :final_in_dim].copy_(stock_head.weight)
            else:
                h.head.weight.copy_(stock_head.weight)
            h.head.bias.copy_(stock_head.bias)

        return h
