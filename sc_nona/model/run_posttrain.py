"""Real-data posttraining entry point for NTv3WithRNA.

In NTv3 terminology this is *posttraining*: we start from an `NTv3_*_pre`
checkpoint (nucleotide-MLM-pretrained) and add new modalities + heads + a
joint objective on top. Same role as NTv3's own posttraining (which adds
bigwig regression heads), just with an RNA *input* stream and masked
reconstruction instead of supervised regression.

Outer-loop structure is borrowed from NTv3's posttraining notebook
(notebooks_tutorials/03_fine_tuning_posttrained_model_biwig.ipynb):
  - step-based (not epoch-based) loop
  - gradient accumulation: K micro-batches per optimizer.step
  - cycling train iterator (auto-restart on StopIteration)
  - periodic train-metric logging + periodic eval-on-val + checkpoint save
  - AdamW + LambdaLR (NTv3's modified-sqrt decay, ported verbatim)

What we don't borrow: NTv3's TracksMetrics (per-track Pearson) and
poisson_multinomial_loss -- those are for bigwig regression. We use plain MSE
on log1p RNA targets + CE on DNA. Per-track validation metrics can be added
later (see TODO in `validate`).

Usage:
    python run_posttrain.py \\
        --fasta /path/to/genome.fasta \\
        --bigwigs cell1.bigwig cell2.bigwig ... \\
        --splits-bed /path/to/splits.bed \\
        --window-length 32768 \\
        --batch-size 4 \\
        --num-steps 100000 \\
        --output-dir checkpoints/

Set HF_TOKEN in the environment so the NTv3 checkpoint can be downloaded.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import torch

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.data import (
    GenomicWindowDataset,
    enumerate_windows,
    make_dataloader,
    read_splits_bed,
)
from model.masking import BlockAlongLengthMasker
from model.model import NTv3WithRNA
from model.train import (
    build_optimizer,
    joint_loss,
    ntv3_sqrt_decay_schedule,
    prepare_batch,
)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


@dataclass
class TrainConfig:
    # Data
    fasta: str
    bigwigs: list[str]
    splits_bed: str
    window_length: int = 32_768
    train_split: str = "train"
    val_split: str = "val"
    stride: int | None = None  # defaults to window_length (non-overlapping tiles)

    # Model
    model_name: str = "InstaDeepAI/NTv3_8M_pre"
    d_rna: int = 64
    num_axial_layers: int = 4
    num_heads: int = 4

    # Optim
    batch_size: int = 4
    grad_accum: int = 1
    num_steps: int = 100_000
    warmup_steps: int = 2_000
    lr: float = 3e-4
    weight_decay: float = 0.01
    grad_clip: float = 1.0
    initial_lr: float = 0.0   # warmup starting LR
    final_ratio: float = 0.5  # LR multiplier at last step

    # Loss weights
    alpha: float = 1.0   # L_dna_stock
    beta: float = 1.0    # L_dna_with_rna
    gamma: float = 1.0   # L_rna

    # Masking
    mask_frac: float = 0.15

    # Loop
    log_every: int = 50
    validate_every: int = 1_000
    output_dir: str = "checkpoints"
    num_workers: int = 2
    seed: int = 0
    local_files_only: bool = False  # skip HF network checks (use cached only)

    # Wandb (all None = disabled; set wandb_project to enable).
    wandb_project: str | None = None
    wandb_entity: str | None = None
    wandb_run_name: str | None = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _cycle(loader):
    """Yield from `loader` forever, re-iterating when exhausted."""
    while True:
        for b in loader:
            yield b


def _format_sample(chrom: str, start: int, end: int,
                   dna_bin_mask_row: torch.Tensor, bin_size: int) -> str:
    """One-line summary of a single training window + the bin range that was masked.

    Works for any masker; for non-contiguous masks (e.g. random-entry) the shown
    bp range is just min(masked_bin)..max(masked_bin)+1 -- coarse but informative.
    """
    masked_idx = dna_bin_mask_row.nonzero(as_tuple=False).flatten()
    win = f"{chrom}:{start:,}-{end:,}"
    if masked_idx.numel() == 0:
        return f"    {win}  mask: (none)"
    bin_lo = int(masked_idx.min())
    bin_hi = int(masked_idx.max()) + 1
    bp_lo = start + bin_lo * bin_size
    bp_hi = start + bin_hi * bin_size
    n_bins = int(dna_bin_mask_row.sum())
    return (
        f"    {win}  "
        f"mask: {chrom}:{bp_lo:,}-{bp_hi:,}  "
        f"({n_bins} bin{'s' if n_bins != 1 else ''}, {bp_hi - bp_lo:,}bp)"
    )


@dataclass
class RunningAvg:
    sums: dict[str, float] = field(default_factory=dict)
    n: int = 0

    def update(self, components: dict[str, float]) -> None:
        for k, v in components.items():
            self.sums[k] = self.sums.get(k, 0.0) + v
        self.n += 1

    def snapshot(self) -> dict[str, float]:
        return {k: v / max(self.n, 1) for k, v in self.sums.items()}

    def reset(self) -> None:
        self.sums = {}
        self.n = 0


def validate(
    model: NTv3WithRNA,
    val_loader,
    masker: BlockAlongLengthMasker,
    *,
    bin_size: int,
    mask_token_id: int,
    device: torch.device,
    cfg: TrainConfig,
) -> dict[str, float]:
    """Run one pass over val_loader; return averaged loss components.

    TODO: per-track Pearson correlation on RNA predictions would be a useful
    additional metric (mirrors NTv3 TracksMetrics).
    """
    model.eval()
    avg = RunningAvg()
    with torch.no_grad():
        for batch in val_loader:
            dna = batch["dna_tokens"].to(device)
            rna = batch["rna_matrix"].to(device)
            step = prepare_batch(
                dna, rna, masker=masker, bin_size=bin_size,
                mask_token_id=mask_token_id, device=device,
            )
            out = model(step.dna_tokens, step.rna_matrix, rna_mask=step.rna_mask)
            _, comps = joint_loss(
                out, step.dna_targets, step.rna_targets, step.rna_mask,
                alpha=cfg.alpha, beta=cfg.beta, gamma=cfg.gamma,
            )
            avg.update(comps)
    model.train()
    return avg.snapshot()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(cfg: TrainConfig) -> None:
    torch.manual_seed(cfg.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"device={device}  output={out_dir}")

    # ---- Wandb (opt-in) ----
    wandb_run = None
    if cfg.wandb_project:
        import wandb
        wandb_run = wandb.init(
            project=cfg.wandb_project,
            entity=cfg.wandb_entity,
            name=cfg.wandb_run_name,
            config=cfg.__dict__,
            dir=str(out_dir),
        )
        print(f"wandb: project={cfg.wandb_project} run={wandb_run.name}")

    # ---- Tokenizer + model ----
    from transformers import AutoTokenizer
    tokenizer = AutoTokenizer.from_pretrained(
        cfg.model_name, trust_remote_code=True,
        token=os.environ.get("HF_TOKEN"),
        local_files_only=cfg.local_files_only,
    )
    model = NTv3WithRNA.from_pretrained(
        cfg.model_name,
        d_rna=cfg.d_rna,
        num_axial_layers=cfg.num_axial_layers,
        num_heads=cfg.num_heads,
        token=os.environ.get("HF_TOKEN"),
        local_files_only=cfg.local_files_only,
    ).to(device)
    bin_size = model.bin_size
    mask_token_id = model.config.mask_token_id
    print(f"model: embed_dim={model.config.embed_dim}, "
          f"num_layers={model.config.num_layers}, bin_size={bin_size}")

    # ---- Data ----
    train_intervals = read_splits_bed(cfg.splits_bed, cfg.train_split)
    val_intervals = read_splits_bed(cfg.splits_bed, cfg.val_split)
    train_windows = enumerate_windows(train_intervals, cfg.window_length, cfg.stride)
    val_windows = enumerate_windows(val_intervals, cfg.window_length, cfg.window_length)
    print(f"train windows: {len(train_windows)}, val windows: {len(val_windows)}")

    train_ds = GenomicWindowDataset(cfg.fasta, cfg.bigwigs, train_windows, tokenizer, bin_size)
    val_ds = GenomicWindowDataset(cfg.fasta, cfg.bigwigs, val_windows, tokenizer, bin_size)
    train_loader = make_dataloader(train_ds, cfg.batch_size, shuffle=True, num_workers=cfg.num_workers)
    val_loader = make_dataloader(val_ds, cfg.batch_size, shuffle=False, num_workers=cfg.num_workers)
    train_iter = _cycle(train_loader)

    # ---- Optim ----
    optimizer = build_optimizer(model, lr=cfg.lr, weight_decay=cfg.weight_decay)
    scheduler = ntv3_sqrt_decay_schedule(
        optimizer, warmup_steps=cfg.warmup_steps, total_steps=cfg.num_steps,
        initial_lr=cfg.initial_lr, peak_lr=cfg.lr, final_ratio=cfg.final_ratio,
    )
    masker = BlockAlongLengthMasker(frac=cfg.mask_frac)

    print(f"starting training for {cfg.num_steps} steps "
          f"(grad_accum={cfg.grad_accum}, effective batch = {cfg.batch_size * cfg.grad_accum})")
    running = RunningAvg()
    best_val_total = math.inf
    t0 = time.time()
    model.train()
    # Track the most recent micro-batch's windows + mask so we can print them at log_every.
    last_meta: dict | None = None

    for step_idx in range(cfg.num_steps):
        optimizer.zero_grad(set_to_none=True)

        # Gradient accumulation: K micro-batches per optimizer.step.
        for _ in range(cfg.grad_accum):
            batch = next(train_iter)
            dna = batch["dna_tokens"].to(device)
            rna = batch["rna_matrix"].to(device)
            step = prepare_batch(
                dna, rna, masker=masker, bin_size=bin_size,
                mask_token_id=mask_token_id, device=device,
            )
            out = model(step.dna_tokens, step.rna_matrix, rna_mask=step.rna_mask)
            loss, comps = joint_loss(
                out, step.dna_targets, step.rna_targets, step.rna_mask,
                alpha=cfg.alpha, beta=cfg.beta, gamma=cfg.gamma,
            )
            (loss / cfg.grad_accum).backward()
            running.update(comps)
            last_meta = {
                "chrom": batch["chrom"],
                "start": batch["start"],
                "end": batch["end"],
                "dna_bin_mask": step.dna_bin_mask.detach().cpu(),
            }

        if cfg.grad_clip is not None:
            torch.nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip)
        optimizer.step()
        scheduler.step()

        # Per-step wandb logging (use the most recent micro-batch's components).
        if wandb_run is not None:
            wandb_run.log(
                {
                    "train/total": comps["total"],
                    "train/l_dna_stock": comps["l_dna_stock"],
                    "train/l_dna_with_rna": comps["l_dna_with_rna"],
                    "train/l_rna": comps["l_rna"],
                    "train/lr": optimizer.param_groups[0]["lr"],
                },
                step=step_idx + 1,
            )

        # Periodic train-loss logging.
        if (step_idx + 1) % cfg.log_every == 0:
            avg = running.snapshot()
            lr_now = optimizer.param_groups[0]["lr"]
            elapsed = time.time() - t0
            print(
                f"step {step_idx + 1:>6}/{cfg.num_steps}  "
                f"lr={lr_now:.2e}  "
                f"total={avg['total']:.4f}  "
                f"dna_stock={avg['l_dna_stock']:.4f}  "
                f"dna_with_rna={avg['l_dna_with_rna']:.4f}  "
                f"rna={avg['l_rna']:.4f}  "
                f"({elapsed:.0f}s)"
            )
            if last_meta is not None:
                for b in range(len(last_meta["chrom"])):
                    print(_format_sample(
                        last_meta["chrom"][b],
                        last_meta["start"][b],
                        last_meta["end"][b],
                        last_meta["dna_bin_mask"][b],
                        bin_size,
                    ))
            running.reset()

        # Periodic validation + checkpoint.
        if (step_idx + 1) % cfg.validate_every == 0:
            print(f"\nvalidating at step {step_idx + 1}...")
            val_metrics = validate(
                model, val_loader, masker,
                bin_size=bin_size, mask_token_id=mask_token_id,
                device=device, cfg=cfg,
            )
            print(
                f"  VAL  total={val_metrics['total']:.4f}  "
                f"dna_stock={val_metrics['l_dna_stock']:.4f}  "
                f"dna_with_rna={val_metrics['l_dna_with_rna']:.4f}  "
                f"rna={val_metrics['l_rna']:.4f}"
            )
            if wandb_run is not None:
                wandb_run.log(
                    {f"val/{k}": v for k, v in val_metrics.items()},
                    step=step_idx + 1,
                )
            last_path = out_dir / f"ckpt_step_{step_idx + 1}.pt"
            torch.save({
                "step": step_idx + 1,
                "model_state": model.state_dict(),
                "optim_state": optimizer.state_dict(),
                "sched_state": scheduler.state_dict(),
                "config": cfg.__dict__,
                "val_metrics": val_metrics,
            }, last_path)
            # Track best-by-val-total; save a separate symlink/copy.
            if val_metrics["total"] < best_val_total:
                best_val_total = val_metrics["total"]
                best_path = out_dir / "best.pt"
                torch.save(torch.load(last_path), best_path)
                print(f"  new best val total {best_val_total:.4f} -> {best_path}")
                if wandb_run is not None:
                    wandb_run.summary["best_val_total"] = best_val_total
                    wandb_run.summary["best_step"] = step_idx + 1
            print()
            model.train()

    if wandb_run is not None:
        wandb_run.finish()
    print("training complete")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args() -> TrainConfig:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--bigwigs", nargs="+", required=True)
    p.add_argument("--splits-bed", required=True)
    p.add_argument("--window-length", type=int, default=32_768)
    p.add_argument("--train-split", default="train")
    p.add_argument("--val-split", default="val")
    p.add_argument("--stride", type=int, default=None)
    p.add_argument("--model-name", default="InstaDeepAI/NTv3_8M_pre")
    p.add_argument("--d-rna", type=int, default=64)
    p.add_argument("--num-axial-layers", type=int, default=4)
    p.add_argument("--num-heads", type=int, default=4)
    p.add_argument("--batch-size", type=int, default=4)
    p.add_argument("--grad-accum", type=int, default=1)
    p.add_argument("--num-steps", type=int, default=100_000)
    p.add_argument("--warmup-steps", type=int, default=2_000)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--weight-decay", type=float, default=0.01)
    p.add_argument("--grad-clip", type=float, default=1.0)
    p.add_argument("--initial-lr", type=float, default=0.0)
    p.add_argument("--final-ratio", type=float, default=0.5)
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--beta", type=float, default=1.0)
    p.add_argument("--gamma", type=float, default=1.0)
    p.add_argument("--mask-frac", type=float, default=0.15)
    p.add_argument("--log-every", type=int, default=50)
    p.add_argument("--validate-every", type=int, default=1_000)
    p.add_argument("--output-dir", default="checkpoints")
    p.add_argument("--num-workers", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--local-files-only", action="store_true",
                   help="Skip HF network checks; use only cached files. "
                        "Use after a successful online run to go fully offline.")
    p.add_argument("--wandb-project", default=None,
                   help="Wandb project. If set, logs train/val losses + LR to wandb.")
    p.add_argument("--wandb-entity", default=None)
    p.add_argument("--wandb-run-name", default=None)
    args = p.parse_args()
    return TrainConfig(
        fasta=args.fasta,
        bigwigs=args.bigwigs,
        splits_bed=args.splits_bed,
        window_length=args.window_length,
        train_split=args.train_split,
        val_split=args.val_split,
        stride=args.stride,
        model_name=args.model_name,
        d_rna=args.d_rna,
        num_axial_layers=args.num_axial_layers,
        num_heads=args.num_heads,
        batch_size=args.batch_size,
        grad_accum=args.grad_accum,
        num_steps=args.num_steps,
        warmup_steps=args.warmup_steps,
        lr=args.lr,
        weight_decay=args.weight_decay,
        grad_clip=args.grad_clip,
        initial_lr=args.initial_lr,
        final_ratio=args.final_ratio,
        alpha=args.alpha,
        beta=args.beta,
        gamma=args.gamma,
        mask_frac=args.mask_frac,
        log_every=args.log_every,
        validate_every=args.validate_every,
        output_dir=args.output_dir,
        num_workers=args.num_workers,
        seed=args.seed,
        local_files_only=args.local_files_only,
        wandb_project=args.wandb_project,
        wandb_entity=args.wandb_entity,
        wandb_run_name=args.wandb_run_name,
    )


if __name__ == "__main__":
    main(_parse_args())
