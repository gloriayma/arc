"""Joint pretraining loop for NTv3WithRNA.

LOSS = alpha * L_dna_stock + beta * L_dna_with_rna + gamma * L_rna_mlm

L_dna_stock     : F.cross_entropy on out["dna_logits_stock"] vs masked dna targets
                  (NTv3 native MLM, conditioned on DNA only)
L_dna_with_rna  : F.cross_entropy on out["dna_logits_with_rna"] vs same targets
                  (DNARNAHead: DNA + RNA -> tokens)
L_rna_mlm       : MSE between out["rna_pred"] and log1p(target) at masked (cell, bin) entries
                  (RNAMlmHead: per-cell continuous reconstruction)
                  TODO: this loss is at bin-level (L/128). Consider switching to
                  1bp resolution later -- mirrors NTv3's bigwig head and would
                  be needed for splice / single-base-resolution downstream tasks.
                  See model/rna_mlm_head.py for the architectural options.

Mask coupling: a `Masker` (model/masking.py) chooses the block(s) once per batch.
For the default BlockAlongLengthMasker, the same L-block is masked in both modalities,
so the RNA-MLM teaches DNA the masked-region DNA->RNA dependency.

Alternative RNA-MLM losses (not implemented; we use MSE-on-log1p for v1):
    * clip + sqrt: clip raw counts at a high quantile (e.g. 99th percentile) and
      sqrt-transform; predict the sqrt'd values with MSE. Variance-stabilizing,
      used by Borzoi. Trades log1p's gentleness toward zeros for less aggressive
      compression at the top of the range.
    * Poisson NLL: softplus the model output to a positive rate `lambda`, loss is
      `lambda - target * log(lambda)`. Natural for raw counts; outlier-sensitive.
    * Poisson-multinomial (Enformer / NTv3 bigwig finetune convention): decompose
      into (a) Poisson on per-track summed counts in the window ("scale") and
      (b) multinomial CE between normalized per-position distributions ("shape").
      Standard weighting is 5x on the shape term. Captures "where the peaks are"
      separately from "how much total signal".
    To switch: change `rna_loss_fn` below and update the RNAEncoder/RNAMlmHead
    preprocessing accordingly (drop log1p, add softplus on outputs, etc.).

NTv3 components we mirror (cf. 06_NTv3_generative_training.ipynb):
    * AdamW over only-trainable params
    * LambdaLR with linear warmup then power-law decay (we use a simpler cosine
      decay; modified_sqrt_decay can be swapped in)
    * Outer loop: zero_grad -> forward -> loss -> backward -> step -> sched.step -> log
    * Gradient accumulation supported by accumulating loss and stepping every K iters

The "overfit one fixed batch" sanity check lives in
model/tests/test_overfit.py.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import torch
from torch import nn
from torch.nn import functional as F
from torch.optim import AdamW
from torch.optim.lr_scheduler import LambdaLR

from model.masking import Masker, apply_dna_mask, rna_targets_from_matrix
from model.model import NTv3WithRNA


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------


def joint_loss(
    out: dict[str, torch.Tensor],
    dna_targets: torch.Tensor,
    rna_targets: torch.Tensor,
    rna_mask: torch.Tensor,
    *,
    alpha: float = 1.0,
    beta: float = 1.0,
    gamma: float = 1.0,
    ignore_index: int = -100,
) -> tuple[torch.Tensor, dict[str, float]]:
    """Compute the three-term joint loss + a dict of scalar components for logging.

    Args:
        out:         NTv3WithRNA.forward output dict
        dna_targets: (B, L) int, original token ids at masked positions, ignore_index elsewhere
        rna_targets: (B, L_bin, D) float, log1p targets at every (cell, bin) entry
        rna_mask:    (B, L_bin, D) bool, True at entries to include in the RNA loss
    """
    # DNA stock: (B, L, alphabet) -> flatten for CE.
    l_dna_stock = F.cross_entropy(
        out["dna_logits_stock"].reshape(-1, out["dna_logits_stock"].shape[-1]),
        dna_targets.reshape(-1),
        ignore_index=ignore_index,
    )
    l_dna_with_rna = F.cross_entropy(
        out["dna_logits_with_rna"].reshape(-1, out["dna_logits_with_rna"].shape[-1]),
        dna_targets.reshape(-1),
        ignore_index=ignore_index,
    )

    # RNA: MSE in log1p space, only at masked entries.
    if rna_mask.any():
        rna_pred_at_mask = out["rna_pred"][rna_mask]
        rna_target_at_mask = rna_targets[rna_mask]
        l_rna = F.mse_loss(rna_pred_at_mask, rna_target_at_mask)
    else:
        l_rna = torch.zeros((), device=out["rna_pred"].device, dtype=out["rna_pred"].dtype)

    total = alpha * l_dna_stock + beta * l_dna_with_rna + gamma * l_rna
    components = {
        "l_dna_stock": float(l_dna_stock.detach()),
        "l_dna_with_rna": float(l_dna_with_rna.detach()),
        "l_rna": float(l_rna.detach()),
        "total": float(total.detach()),
    }
    return total, components


# ---------------------------------------------------------------------------
# Optimizer + scheduler
# ---------------------------------------------------------------------------


def build_optimizer(
    model: nn.Module,
    lr: float = 3e-4,
    weight_decay: float = 0.01,
    betas: tuple[float, float] = (0.9, 0.95),
) -> AdamW:
    params = [p for p in model.parameters() if p.requires_grad]
    return AdamW(params, lr=lr, betas=betas, weight_decay=weight_decay)


def cosine_warmup_schedule(
    optimizer: torch.optim.Optimizer,
    warmup_steps: int,
    total_steps: int,
    min_lr_ratio: float = 0.1,
) -> LambdaLR:
    """Linear warmup then cosine decay to `min_lr_ratio * peak_lr`.

    (NTv3 uses linear warmup + power-law/sqrt decay -- see
    `get_modified_sqrt_decay_scheduler` in 06_NTv3_generative_training.ipynb;
    cosine is a fine substitute for our purposes and easier to reason about.)
    """
    def lr_lambda(step: int) -> float:
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        cosine = 0.5 * (1.0 + math.cos(math.pi * min(progress, 1.0)))
        return min_lr_ratio + (1.0 - min_lr_ratio) * cosine

    return LambdaLR(optimizer, lr_lambda=lr_lambda)


# ---------------------------------------------------------------------------
# Step
# ---------------------------------------------------------------------------


@dataclass
class StepInputs:
    dna_tokens: torch.Tensor      # (B, L) masked input ids (with mask_token_id at masked bins)
    rna_matrix: torch.Tensor      # (B, L, D) raw RNA coverage (un-masked input for the encoder)
    dna_targets: torch.Tensor     # (B, L) original token ids at masked positions, ignore_index elsewhere
    rna_targets: torch.Tensor     # (B, L_bin, D) log1p targets
    rna_mask: torch.Tensor        # (B, L_bin, D) bool, True at masked entries


def prepare_batch(
    raw_dna_tokens: torch.Tensor,   # (B, L) unmasked
    raw_rna_matrix: torch.Tensor,   # (B, L, D)
    *,
    masker: Masker,
    bin_size: int,
    mask_token_id: int,
    device: torch.device,
    generator: torch.Generator | None = None,
) -> StepInputs:
    """Apply masking and produce all tensors the train_step needs."""
    B, L = raw_dna_tokens.shape
    _, _, D = raw_rna_matrix.shape
    L_bin = L // bin_size

    dna_bin_mask, rna_mask = masker(B, L_bin, D, device=device, generator=generator)
    masked_input_ids, dna_targets = apply_dna_mask(
        raw_dna_tokens, dna_bin_mask, bin_size=bin_size, mask_token_id=mask_token_id
    )
    rna_targets = rna_targets_from_matrix(raw_rna_matrix, bin_size)
    return StepInputs(
        dna_tokens=masked_input_ids,
        rna_matrix=raw_rna_matrix,
        dna_targets=dna_targets,
        rna_targets=rna_targets,
        rna_mask=rna_mask,
    )


def train_step(
    model: NTv3WithRNA,
    batch: StepInputs,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler._LRScheduler | None,
    *,
    alpha: float = 1.0,
    beta: float = 1.0,
    gamma: float = 1.0,
    grad_clip: float | None = 1.0,
) -> dict[str, float]:
    model.train()
    optimizer.zero_grad(set_to_none=True)

    out = model(batch.dna_tokens, batch.rna_matrix, rna_mask=batch.rna_mask)
    loss, components = joint_loss(
        out, batch.dna_targets, batch.rna_targets, batch.rna_mask,
        alpha=alpha, beta=beta, gamma=gamma,
    )
    loss.backward()
    if grad_clip is not None:
        torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)
    optimizer.step()
    if scheduler is not None:
        scheduler.step()
    components["lr"] = optimizer.param_groups[0]["lr"]
    return components


