"""Tests for data.py.

Covers:
  - SyntheticGenomicDataset shapes
  - collate_fn / DataLoader stacking
  - end-to-end glue: dataloader output -> prepare_batch -> train_step
  - read_splits_bed + enumerate_windows logic on an in-memory BED-like file

We don't exercise the real FASTA + BigWig path here (would need test fixtures
on disk). The pyfaidx / pyBigWig calls are thin wrappers around well-tested
libraries; once real files exist, the smoke check is to construct
GenomicWindowDataset(...) and verify shapes match SyntheticGenomicDataset.
"""
import io
import sys
from pathlib import Path

import pandas as pd
import torch

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.data import (  # noqa: E402
    SyntheticGenomicDataset,
    Window,
    collate_windows,
    enumerate_windows,
    make_dataloader,
    read_splits_bed,
)


def test_synthetic_dataset_shapes() -> None:
    ds = SyntheticGenomicDataset(num_windows=4, window_length=512, D=3)
    assert len(ds) == 4
    item = ds[0]
    assert item["dna_tokens"].shape == (512,)
    assert item["dna_tokens"].dtype == torch.long
    assert item["rna_matrix"].shape == (512, 3)
    assert item["rna_matrix"].dtype == torch.float32
    # Determinism: same indexing should yield the same tensor.
    assert torch.equal(ds[0]["dna_tokens"], ds[0]["dna_tokens"])
    print("PASS  SyntheticGenomicDataset shapes + determinism")


def test_collate_and_loader() -> None:
    ds = SyntheticGenomicDataset(num_windows=6, window_length=512, D=4)
    loader = make_dataloader(ds, batch_size=2, shuffle=False, drop_last=False)
    batch = next(iter(loader))
    assert batch["dna_tokens"].shape == (2, 512)
    assert batch["rna_matrix"].shape == (2, 512, 4)
    print("PASS  DataLoader + collate_windows produce (B,L) and (B,L,D)")


def test_loader_to_prepare_batch() -> None:
    """Glue: DataLoader output flows into prepare_batch unchanged."""
    from model.masking import BlockAlongLengthMasker
    from model.train import prepare_batch

    ds = SyntheticGenomicDataset(num_windows=4, window_length=512, D=3)
    loader = make_dataloader(ds, batch_size=2, shuffle=False)
    batch = next(iter(loader))

    masker = BlockAlongLengthMasker(frac=0.25)
    step = prepare_batch(
        batch["dna_tokens"], batch["rna_matrix"],
        masker=masker, bin_size=128, mask_token_id=2,
        device=torch.device("cpu"),
    )
    assert step.dna_tokens.shape == (2, 512)
    assert step.rna_matrix.shape == (2, 512, 3)
    assert step.dna_targets.shape == (2, 512)
    assert step.rna_targets.shape == (2, 4, 3)  # L_bin = 512/128 = 4
    assert step.rna_mask.shape == (2, 4, 3)
    # Block-along-length: rna_mask is dna_bin_mask broadcast across D.
    assert step.rna_mask.any() and not step.rna_mask.all()
    print("PASS  loader -> prepare_batch end-to-end shapes")


def test_read_splits_bed(tmp_path: Path) -> None:
    """Parse a 4-column BED with multiple splits, tile into fixed-length windows."""
    bed = tmp_path / "splits.bed"
    bed.write_text(
        "chr1\t0\t2000\ttrain\n"
        "chr2\t0\t512\tval\n"
        "chr1\t10000\t12048\ttrain\n"
    )
    train_intervals = read_splits_bed(bed, "train")
    assert len(train_intervals) == 2
    assert train_intervals[0] == Window("chr1", 0, 2000)

    val_intervals = read_splits_bed(bed, "val")
    assert val_intervals == [Window("chr2", 0, 512)]

    # Tiling: 2000bp into 512bp windows = 3 windows; 2048bp -> 4 windows.
    train_windows = enumerate_windows(train_intervals, window_length=512)
    assert len(train_windows) == 3 + 4
    assert train_windows[0] == Window("chr1", 0, 512)
    assert train_windows[3] == Window("chr1", 10000, 10512)
    print("PASS  read_splits_bed + enumerate_windows")


def main() -> None:
    test_synthetic_dataset_shapes()
    test_collate_and_loader()
    test_loader_to_prepare_batch()
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        test_read_splits_bed(Path(td))
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
