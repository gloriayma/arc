"""Overfit-one-batch sanity check for the joint train loop.

Fixes a single synthetic batch and a single mask, then trains repeatedly,
asserting that each of the three loss components (DNA-stock, DNA-with-RNA,
RNA-MLM) drops by at least 2x from start to end. This exercises:

  - the full NTv3WithRNA forward
  - all three loss terms in joint_loss
  - gradient flow into core, RNAEncoder, RNAMlmHead, DNARNAHead
  - the optimizer + scheduler wiring in train.py

Run with HF_TOKEN set so the NTv3 checkpoint can be downloaded.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import torch

_REPO = Path(__file__).resolve().parents[2]  # sc_nona/
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.masking import BlockAlongLengthMasker  # noqa: E402
from model.model import NTv3WithRNA  # noqa: E402
from model.train import (  # noqa: E402
    build_optimizer,
    ntv3_sqrt_decay_schedule,
    prepare_batch,
    train_step,
)


def _make_synthetic_batch(
    B: int, L: int, D: int, alphabet_size: int, device: torch.device, seed: int = 0,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Random DNA tokens in [3, alphabet_size) (skip pad=1, mask=2) and random RNA."""
    g = torch.Generator(device=device).manual_seed(seed)
    dna = torch.randint(low=3, high=alphabet_size, size=(B, L), device=device, generator=g)
    rna = torch.rand(B, L, D, device=device, generator=g) * 10
    return dna, rna


def overfit_one_batch(
    model_name: str = "InstaDeepAI/NTv3_8M_pre",
    *,
    B: int = 2, L: int = 512, D: int = 4,
    d_rna: int = 32, num_axial_layers: int = 2,
    steps: int = 200, log_every: int = 25,
    lr: float = 3e-4,
    token: str | None = None,
) -> dict[str, list[float]]:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}")

    model = NTv3WithRNA.from_pretrained(
        model_name, d_rna=d_rna, num_axial_layers=num_axial_layers, token=token,
    ).to(device)
    cfg = model.config
    mask_token_id = cfg.mask_token_id
    bin_size = model.bin_size

    raw_dna, raw_rna = _make_synthetic_batch(B, L, D, cfg.alphabet_size, device)
    masker = BlockAlongLengthMasker(frac=0.25)
    # Fix the mask too so we're truly overfitting one batch.
    g = torch.Generator(device=device).manual_seed(123)
    batch = prepare_batch(
        raw_dna, raw_rna, masker=masker,
        bin_size=bin_size, mask_token_id=mask_token_id,
        device=device, generator=g,
    )

    optimizer = build_optimizer(model, lr=lr)
    scheduler = ntv3_sqrt_decay_schedule(
        optimizer, warmup_steps=20, total_steps=steps, peak_lr=lr,
    )

    history: dict[str, list[float]] = {
        "l_dna_stock": [], "l_dna_with_rna": [], "l_rna": [], "total": [], "lr": []
    }
    for step in range(steps):
        comps = train_step(
            model, batch, optimizer, scheduler,
            alpha=1.0, beta=1.0, gamma=1.0,
        )
        for k, v in comps.items():
            history[k].append(v)
        if step % log_every == 0 or step == steps - 1:
            print(
                f"step {step:4d}  total={comps['total']:.4f}  "
                f"dna_stock={comps['l_dna_stock']:.4f}  "
                f"dna_with_rna={comps['l_dna_with_rna']:.4f}  "
                f"rna={comps['l_rna']:.4f}  "
                f"lr={comps['lr']:.2e}"
            )

    for k in ("l_dna_stock", "l_dna_with_rna", "l_rna"):
        first = sum(history[k][:5]) / 5
        last = sum(history[k][-5:]) / 5
        print(f"  {k}: avg first 5 = {first:.4f}, avg last 5 = {last:.4f}")
        assert last < first * 0.5, f"{k} did not decrease by 2x: {first:.4f} -> {last:.4f}"

    print("\nOVERFIT-ONE-BATCH: PASS")
    return history


if __name__ == "__main__":
    token = os.environ.get("HF_TOKEN")
    overfit_one_batch(token=token)
