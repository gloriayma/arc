"""Masking strategies for joint DNA + scRNA MLM training.

Design: a `Masker` is any callable
    masker(B, L_bin, D, *, device, generator=None) -> (dna_bin_mask, rna_mask)

that returns two boolean tensors:
    dna_bin_mask : (B, L_bin)      True at L/bin bins whose DNA tokens are to be masked
    rna_mask     : (B, L_bin, D)   True at (cell, bin) entries to be predicted

The two masks are returned together so strategies that couple them (e.g. block-
along-length, where the same region is masked in both modalities) can express
that coupling explicitly; strategies that mask only one modality return all-False
for the other.

Helpers below convert these bin-level masks into the per-token / per-loss
tensors the training loop needs:
    apply_dna_mask(input_ids, dna_bin_mask, ...) -> (masked_input_ids, dna_targets)
    rna_targets_from_mask(rna_matrix, rna_mask, bin_size) -> bin-resolution targets

Swapping a strategy = instantiate a different Masker. The train loop is
strategy-agnostic.

What each masker actually teaches:
    BlockAlongLengthMasker masks a contiguous block of bins across *every* cell
    simultaneously, so the model never has to do "predict cell 5's expression at
    bin 100 given that cells 1-4's expression at bin 100 is visible." It teaches
    long-range DNA -> RNA reasoning (from flanking DNA + cell identity in unmasked
    regions), but does NOT directly train cross-cell covariation -- even though
    the axial-D attention in RNAEncoder is structurally capable of learning it.
    A natural next masker is one that drops random (cell, bin) entries independently
    (or whole cells but not whole bins): that's the masker that exploits the D-axis
    attention and teaches cell-cell covariation. We probably want both, mixed per
    batch, once we get past v1.

    Block masking vs. standard BERT scatter-masking: scatter-masking lets the model
    interpolate from local context, which is the easier task; block masking forces
    long-range reasoning since there's no local context inside the masked region.

TODO: think about center-crop symmetry. BlockAlongLengthMasker currently lets
    `starts = randint(0, L_bin - width)`, so the mask can land at the very edge.
    No conv-receptive-field issue on the RNA side (pure attention), but the DNA
    core has conv stems -- edge tokens see asymmetric context. If long-range
    flanking signal turns out to matter, clamp `max_start` to leave some margin
    on each side (NTv3-bigwig style: supervise the middle, give context on both
    flanks). Probably overkill for v1 but worth revisiting.
"""
from __future__ import annotations

from typing import Protocol

import torch
from torch.nn import functional as F


# ---------------------------------------------------------------------------
# Masker protocol + concrete strategies
# ---------------------------------------------------------------------------


class Masker(Protocol):
    def __call__(
        self,
        B: int,
        L_bin: int,
        D: int,
        *,
        device: torch.device,
        generator: torch.Generator | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]: ...


class BlockAlongLengthMasker:
    """Mask one contiguous L-block per example, applied to both DNA and RNA.

    For each batch element, sample a contiguous range of L/bin bins (random start,
    fixed width = round(frac * L_bin), at least min_bins). The same range is masked
    in DNA (all tokens) and RNA (all cells x masked bins).

    Use this when you want to teach the model long-range DNA -> RNA reasoning:
    in the masked region the model has only flanking DNA + cell identity from
    unmasked regions to work from.

    TODO: also try masking only on the RNA. or only on some RNA tracks. 
    """

    def __init__(self, frac: float = 0.15, min_bins: int = 1):
        assert 0.0 < frac < 1.0
        self.frac = frac
        self.min_bins = min_bins

    def __call__(
        self,
        B: int,
        L_bin: int,
        D: int,
        *,
        device: torch.device,
        generator: torch.Generator | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        width = max(self.min_bins, round(self.frac * L_bin))
        width = min(width, L_bin)
        max_start = L_bin - width
        # Random start per example.
        starts = torch.randint(
            low=0, high=max_start + 1, size=(B,),
            device=device, generator=generator,
        )

        # Build the (B, L_bin) mask via broadcasting.
        idx = torch.arange(L_bin, device=device).unsqueeze(0)        # (1, L_bin)
        dna_bin_mask = (idx >= starts.unsqueeze(1)) & (idx < (starts + width).unsqueeze(1))

        # RNA mask = DNA mask broadcast across all D cells.
        rna_mask = dna_bin_mask.unsqueeze(-1).expand(B, L_bin, D).contiguous()
        return dna_bin_mask, rna_mask


# ---------------------------------------------------------------------------
# Helpers to apply masks at training time
# ---------------------------------------------------------------------------


def apply_dna_mask(
    input_ids: torch.Tensor,
    dna_bin_mask: torch.Tensor,
    bin_size: int,
    mask_token_id: int,
    ignore_index: int = -100,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Expand bin-level mask to token-level, then construct masked inputs + CE targets.

    Args:
        input_ids:    (B, L) int token ids
        dna_bin_mask: (B, L_bin) bool from a Masker
        bin_size:     L / L_bin (NTv3: 2 ** num_downsamples = 128)
        mask_token_id: tokenizer mask id (NTv3 pretrained: 2)
        ignore_index: PyTorch CE ignore_index (default -100)

    Returns:
        masked_input_ids: (B, L) input_ids with mask_token_id at masked positions
        dna_targets:      (B, L) original tokens at masked positions, ignore_index elsewhere
    """
    B, L = input_ids.shape
    _, L_bin = dna_bin_mask.shape
    assert L == L_bin * bin_size, f"L={L} != L_bin={L_bin} * bin_size={bin_size}"

    # (B, L_bin) -> (B, L) by repeat_interleave.
    token_mask = dna_bin_mask.repeat_interleave(bin_size, dim=1)  # (B, L)

    masked_input_ids = torch.where(
        token_mask, torch.full_like(input_ids, mask_token_id), input_ids
    )
    dna_targets = torch.where(
        token_mask, input_ids, torch.full_like(input_ids, ignore_index)
    )
    return masked_input_ids, dna_targets


def rna_targets_from_matrix(
    rna_matrix: torch.Tensor,
    bin_size: int,
) -> torch.Tensor:
    """Mean-pool the raw RNA matrix to L/bin and return log1p (the encoder's
    input preprocessing). This is the regression target for the RNA-MLM head,
    which predicts in log1p-space.

    Args:
        rna_matrix: (B, L, D) raw per-base coverage
    Returns:
        (B, L_bin, D) log1p of mean-pooled coverage
    """
    binned = rna_matrix.transpose(1, 2)  # (B, D, L)
    binned = F.avg_pool1d(binned, kernel_size=bin_size)  # (B, D, L_bin)
    binned = binned.transpose(1, 2)  # (B, L_bin, D)
    return torch.log1p(binned.clamp(min=0))
