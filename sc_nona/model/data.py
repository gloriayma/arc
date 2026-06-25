"""Dataset / DataLoader for joint DNA + scRNA training.

Two flavors:

  GenomicWindowDataset       -- real data: FASTA + a list of scRNA BigWigs
                                (one BigWig = one cell/track = one column of the
                                D axis). Window coordinates come from a list of
                                (chrom, start, end) tuples; usually parsed from
                                a `splits.bed` file (NTv3 convention).

  SyntheticGenomicDataset    -- in-memory random tensors that match the real
                                dataset's __getitem__ signature. Useful for
                                pipeline tests without needing real .fasta /
                                .bigwig files on disk.

Both yield dicts:
    {
      "dna_tokens": (L,) int64 token ids,
      "rna_matrix": (L, D) float32 coverage values,
    }

A collate_fn stacks these into (B, L) and (B, L, D) batches; the resulting
DataLoader output can be fed directly to `train.prepare_batch`.

Notes:
  * NaN entries in BigWigs (positions with no read coverage) are replaced with 0.
  * Window lengths must be a multiple of the model's bin_size (NTv3 default 128).
  * No reverse-complement augmentation in v1; add later if useful.
  * No strand handling; scRNA BigWigs are usually un-stranded.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset


# ---------------------------------------------------------------------------
# Window list helpers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Window:
    chrom: str
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive


def read_splits_bed(
    splits_path: str | Path,
    split: str,
    *,
    chrom_col: str = "chr_name",
) -> list[Window]:
    """Parse a 4-column BED-style splits file (chr, start, end, split).

    Matches the NTv3 benchmark convention -- see
    notebooks_tutorials/03_fine_tuning_posttrained_model_biwig.ipynb,
    prepare_genomics_inputs.
    """
    df = pd.read_csv(
        splits_path, sep="\t", header=None,
        names=[chrom_col, "start", "end", "split"],
        dtype={chrom_col: str, "start": int, "end": int, "split": str},
    )
    df = df[df["split"] == split].reset_index(drop=True)
    return [Window(r[chrom_col], int(r["start"]), int(r["end"])) for _, r in df.iterrows()]


def enumerate_windows(
    intervals: Iterable[Window],
    window_length: int,
    stride: int | None = None,
) -> list[Window]:
    """Tile each interval into fixed-length sub-windows of `window_length`.

    stride defaults to window_length (non-overlapping). Sub-windows past the
    interval end are dropped.
    """
    stride = stride or window_length
    out: list[Window] = []
    for iv in intervals:
        s = iv.start
        while s + window_length <= iv.end:
            out.append(Window(iv.chrom, s, s + window_length))
            s += stride
    return out


# ---------------------------------------------------------------------------
# Real dataset (FASTA + BigWigs)
# ---------------------------------------------------------------------------


class GenomicWindowDataset(Dataset):
    """Yields {"dna_tokens", "rna_matrix"} dicts for a list of genomic windows.

    Args:
        fasta_path:     path to a .fasta (indexed via pyfaidx; .fai created on first open)
        bigwig_paths:   list of paths to D scRNA BigWigs (D = cells/tracks)
        windows:        list of Window objects (chrom, start, end) all of length L
        tokenizer:      NTv3 tokenizer (AutoTokenizer.from_pretrained(...))
        bin_size:       must divide L; we assert this once in __init__
    """

    def __init__(
        self,
        fasta_path: str | Path,
        bigwig_paths: Sequence[str | Path],
        windows: Sequence[Window],
        tokenizer,
        bin_size: int = 128,
    ):
        import pyfaidx  # noqa: F401 -- lazy import so the synthetic dataset doesn't need it
        import pyBigWig  # noqa: F401

        self.fasta_path = str(fasta_path)
        self.bigwig_paths = [str(p) for p in bigwig_paths]
        self.windows = list(windows)
        self.tokenizer = tokenizer
        self.bin_size = bin_size
        self.D = len(self.bigwig_paths)
        if not self.windows:
            raise ValueError("no windows provided")
        L = self.windows[0].end - self.windows[0].start
        assert L % bin_size == 0, f"window length {L} not divisible by bin_size {bin_size}"
        for w in self.windows:
            assert (w.end - w.start) == L, f"non-uniform window length: {w}"
        self.window_length = L

        # File handles are lazily opened per worker (multi-process DataLoader safe).
        self._fasta = None
        self._bigwigs: list = []

    def _ensure_open(self) -> None:
        if self._fasta is None:
            import pyfaidx
            self._fasta = pyfaidx.Fasta(self.fasta_path, as_raw=True, sequence_always_upper=True)
        if not self._bigwigs:
            import pyBigWig
            self._bigwigs = [pyBigWig.open(p) for p in self.bigwig_paths]

    def __len__(self) -> int:
        return len(self.windows)

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        self._ensure_open()
        w = self.windows[idx]

        # DNA: read sequence, tokenize.
        seq = str(self._fasta[w.chrom][w.start:w.end])
        # NTv3 tokenizer expects multiple-of-128 length, no special tokens.
        out = self.tokenizer(
            seq, add_special_tokens=False, padding=False, return_tensors="pt",
        )
        dna_tokens = out["input_ids"].squeeze(0).to(torch.long)  # (L,)
        assert dna_tokens.shape[0] == self.window_length, (
            f"tokenized length {dna_tokens.shape[0]} != window_length {self.window_length}"
        )

        # RNA: D BigWigs over the same coords. NaNs -> 0.
        L = self.window_length
        rna = np.zeros((L, self.D), dtype=np.float32)
        for d, bw in enumerate(self._bigwigs):
            vals = bw.values(w.chrom, w.start, w.end, numpy=True)
            vals = np.nan_to_num(vals, nan=0.0, posinf=0.0, neginf=0.0)
            rna[:, d] = vals
        rna_matrix = torch.from_numpy(rna)

        return {"dna_tokens": dna_tokens, "rna_matrix": rna_matrix}


# ---------------------------------------------------------------------------
# Synthetic dataset (no real files needed)
# ---------------------------------------------------------------------------


class SyntheticGenomicDataset(Dataset):
    """Random tensors matching the real dataset signature. For tests / smoke runs."""

    def __init__(
        self,
        num_windows: int,
        window_length: int,
        D: int,
        alphabet_size: int = 11,
        seed: int = 0,
        rna_scale: float = 10.0,
    ):
        self.num_windows = num_windows
        self.window_length = window_length
        self.D = D
        self.alphabet_size = alphabet_size
        self.rna_scale = rna_scale
        # Pre-generate so every epoch sees the same dataset (deterministic).
        g = torch.Generator().manual_seed(seed)
        self._dna = torch.randint(
            low=3, high=alphabet_size, size=(num_windows, window_length),
            dtype=torch.long, generator=g,
        )
        self._rna = (torch.rand(num_windows, window_length, D, generator=g) * rna_scale).float()

    def __len__(self) -> int:
        return self.num_windows

    def __getitem__(self, idx: int) -> dict[str, torch.Tensor]:
        return {"dna_tokens": self._dna[idx], "rna_matrix": self._rna[idx]}


# ---------------------------------------------------------------------------
# Collate + DataLoader factory
# ---------------------------------------------------------------------------


def collate_windows(samples: list[dict[str, torch.Tensor]]) -> dict[str, torch.Tensor]:
    dna_tokens = torch.stack([s["dna_tokens"] for s in samples], dim=0)  # (B, L)
    rna_matrix = torch.stack([s["rna_matrix"] for s in samples], dim=0)  # (B, L, D)
    return {"dna_tokens": dna_tokens, "rna_matrix": rna_matrix}


def make_dataloader(
    dataset: Dataset,
    batch_size: int,
    *,
    shuffle: bool = True,
    num_workers: int = 0,
    drop_last: bool = True,
) -> DataLoader:
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        collate_fn=collate_windows,
        drop_last=drop_last,
        pin_memory=torch.cuda.is_available(),
    )
