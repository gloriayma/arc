"""Tests for the masking module."""
import sys
from pathlib import Path

import torch

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.masking import (  # noqa: E402
    BlockAlongLengthMasker,
    apply_dna_mask,
    rna_targets_from_matrix,
)


def test_block_along_length() -> None:
    torch.manual_seed(0)
    B, L_bin, D = 4, 32, 5
    masker = BlockAlongLengthMasker(frac=0.25, min_bins=1)

    dna_bin_mask, rna_mask = masker(B, L_bin, D, device=torch.device("cpu"))
    assert dna_bin_mask.shape == (B, L_bin)
    assert rna_mask.shape == (B, L_bin, D)
    assert dna_bin_mask.dtype == torch.bool
    assert rna_mask.dtype == torch.bool

    # Each row should have exactly width=round(0.25*32)=8 masked bins.
    expected_width = round(0.25 * L_bin)
    per_row = dna_bin_mask.sum(dim=1)
    assert (per_row == expected_width).all(), f"per_row counts = {per_row.tolist()}"

    # Each row's mask should be a CONTIGUOUS block.
    for b in range(B):
        row = dna_bin_mask[b]
        ones = row.nonzero(as_tuple=False).flatten().tolist()
        assert ones == list(range(ones[0], ones[0] + expected_width)), (
            f"row {b} mask not contiguous: {ones}"
        )

    # RNA mask = DNA mask broadcast across D.
    expected_rna = dna_bin_mask.unsqueeze(-1).expand(B, L_bin, D)
    assert torch.equal(rna_mask, expected_rna), "rna_mask != broadcast(dna_bin_mask)"

    # Different calls should usually pick different start positions.
    dna2, _ = masker(B, L_bin, D, device=torch.device("cpu"))
    assert not torch.equal(dna_bin_mask, dna2), "consecutive calls produced identical mask"
    print("PASS  BlockAlongLengthMasker: shapes, contiguity, broadcast coupling, randomness")


def test_apply_dna_mask() -> None:
    torch.manual_seed(1)
    B, L_bin = 2, 4
    bin_size = 128
    L = L_bin * bin_size
    mask_token_id = 2
    ignore_index = -100

    input_ids = torch.randint(low=3, high=11, size=(B, L), dtype=torch.long)
    # Mask bins 1 and 2 in both rows.
    dna_bin_mask = torch.zeros(B, L_bin, dtype=torch.bool)
    dna_bin_mask[:, 1:3] = True

    masked, targets = apply_dna_mask(input_ids, dna_bin_mask, bin_size, mask_token_id, ignore_index)
    assert masked.shape == input_ids.shape
    assert targets.shape == input_ids.shape

    # Token-level mask: 2 bins x 128 bp = 256 masked tokens per row.
    token_mask = dna_bin_mask.repeat_interleave(bin_size, dim=1)
    assert token_mask.sum(dim=1).tolist() == [256, 256]

    # Masked positions in `masked` are mask_token_id; unmasked retain original.
    assert (masked[token_mask] == mask_token_id).all()
    assert (masked[~token_mask] == input_ids[~token_mask]).all()

    # Targets: original tokens at masked positions, ignore_index elsewhere.
    assert (targets[token_mask] == input_ids[token_mask]).all()
    assert (targets[~token_mask] == ignore_index).all()
    print("PASS  apply_dna_mask: shapes, token-level expansion, target labelling")


def test_rna_targets() -> None:
    B, L, D = 1, 512, 3
    bin_size = 128
    rna = torch.zeros(B, L, D)
    rna[0, :256, 0] = 1.0  # nonzero coverage in first half for cell 0
    targets = rna_targets_from_matrix(rna, bin_size)
    L_bin = L // bin_size
    assert targets.shape == (B, L_bin, D)
    # Bins 0,1 should have nonzero log1p(1)=log(2); bins 2,3 should be 0.
    assert targets[0, 0, 0].item() > 0
    assert targets[0, 1, 0].item() > 0
    assert targets[0, 2, 0].item() == 0
    assert targets[0, 3, 0].item() == 0
    assert (targets[0, :, 1:] == 0).all()  # other cells all zero
    print("PASS  rna_targets_from_matrix: shape + log1p-space binning")


def test_swap_interface() -> None:
    """Demonstrate that swapping a masker is trivial: any object with the
    `masker(B, L_bin, D, *, device, generator=None) -> (dna_bin_mask, rna_mask)`
    signature works."""

    class AllUnmasked:
        def __call__(self, B, L_bin, D, *, device, generator=None):
            return (
                torch.zeros(B, L_bin, dtype=torch.bool, device=device),
                torch.zeros(B, L_bin, D, dtype=torch.bool, device=device),
            )

    for masker in [BlockAlongLengthMasker(frac=0.1), AllUnmasked()]:
        d, r = masker(2, 16, 4, device=torch.device("cpu"))
        assert d.shape == (2, 16)
        assert r.shape == (2, 16, 4)
    print("PASS  Masker protocol: any callable with the right signature is a drop-in")


def main() -> None:
    test_block_along_length()
    test_apply_dna_mask()
    test_rna_targets()
    test_swap_interface()
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
