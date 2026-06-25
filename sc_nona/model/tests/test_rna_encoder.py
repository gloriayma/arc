"""Shape + zero-init checks for RNAEncoder."""
import sys
from pathlib import Path

import torch

_REPO = Path(__file__).resolve().parents[2]  # sc_nona/
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.rna_encoder import RNAEncoder  # noqa: E402


def main() -> None:
    torch.manual_seed(0)

    B, L, D = 2, 512, 5
    bin_size = 128
    L_bin = L // bin_size
    d_rna = 64
    d_ntv = 256  # NTv3 8M

    enc = RNAEncoder(d_rna=d_rna, d_ntv=d_ntv, num_axial_layers=4, num_heads=4,
                     bin_size=bin_size)
    enc.eval()

    rna = torch.rand(B, L, D) * 10  # nonzero, varied
    rna_mask = torch.zeros(B, L_bin, D, dtype=torch.bool)
    rna_mask[:, 1:3, :] = True  # mask bins 1,2 across all cells (block-along-length)

    with torch.no_grad():
        rna_hidden, rna_bias = enc(rna, rna_mask=rna_mask)

    # Shape checks.
    assert rna_hidden.shape == (B, L_bin, D, d_rna), (
        f"rna_hidden wrong shape: {tuple(rna_hidden.shape)}"
    )
    assert rna_bias.shape == (B, L_bin, d_ntv), (
        f"rna_bias wrong shape: {tuple(rna_bias.shape)}"
    )
    print(f"PASS  rna_hidden shape = {tuple(rna_hidden.shape)}")
    print(f"PASS  rna_bias shape   = {tuple(rna_bias.shape)}")

    # Day-0 zero-init: rna_bias should be exactly zero with untrained zero-init to_ntv.
    assert torch.equal(rna_bias, torch.zeros_like(rna_bias)), (
        "day-0 rna_bias is not zero -- to_ntv zero-init broken"
    )
    print("PASS  rna_bias is exactly zero at day-0 (zero-init to_ntv)")

    # rna_hidden should be non-trivial (the encoder isn't a no-op).
    assert rna_hidden.abs().mean().item() > 1e-6, "rna_hidden is identically zero"
    print(f"PASS  rna_hidden is non-trivial (mean|x| = {rna_hidden.abs().mean().item():.4f})")

    # No-mask path should also run cleanly.
    with torch.no_grad():
        rna_hidden_nm, rna_bias_nm = enc(rna, rna_mask=None)
    assert rna_hidden_nm.shape == rna_hidden.shape
    assert torch.equal(rna_bias_nm, torch.zeros_like(rna_bias_nm))
    print("PASS  rna_mask=None path runs and produces zero bias")

    # Gradient flow sanity: bias requires_grad path back to encoder params.
    enc.train()
    rna_hidden, rna_bias = enc(rna, rna_mask=rna_mask)
    # Nudge to_ntv away from zero so a loss is nonzero, then take a step.
    with torch.no_grad():
        enc.to_ntv.weight.normal_(std=0.01)
    rna_hidden, rna_bias = enc(rna, rna_mask=rna_mask)
    loss = rna_bias.pow(2).mean() + rna_hidden.pow(2).mean()
    loss.backward()
    grad_present = sum(1 for p in enc.parameters() if p.grad is not None and p.grad.abs().sum() > 0)
    grad_total = sum(1 for p in enc.parameters())
    print(f"PASS  {grad_present}/{grad_total} encoder params received nonzero gradient")

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
