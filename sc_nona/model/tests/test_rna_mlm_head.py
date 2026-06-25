"""Shape + permutation-equivariance + gradient checks for RNAMlmHead."""
import sys
from pathlib import Path

import torch

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.rna_mlm_head import RNAMlmHead  # noqa: E402


def main() -> None:
    torch.manual_seed(0)

    B, L_bin, D = 2, 4, 5
    d_rna, d_ntv = 64, 256

    head = RNAMlmHead(d_rna=d_rna, d_ntv=d_ntv)
    head.eval()

    rna_hidden = torch.randn(B, L_bin, D, d_rna)
    dna_repr = torch.randn(B, L_bin, d_ntv)

    with torch.no_grad():
        pred = head(rna_hidden, dna_repr)

    # Shape check.
    assert pred.shape == (B, L_bin, D), f"wrong shape {tuple(pred.shape)}"
    print(f"PASS  pred shape = {tuple(pred.shape)}")

    # Permutation equivariance: permuting D in rna_hidden permutes pred in D.
    perm = torch.tensor([2, 0, 4, 1, 3])
    rna_perm = rna_hidden[:, :, perm, :]
    with torch.no_grad():
        pred_perm = head(rna_perm, dna_repr)
    expected = pred[:, :, perm]
    assert torch.allclose(pred_perm, expected, atol=1e-6), (
        "head is not permutation-equivariant along D"
    )
    print("PASS  permutation-equivariant along D (same MLP applied per cell)")

    # Sanity: dna_repr-only path matters (changes in DNA change predictions).
    dna_alt = torch.randn(B, L_bin, d_ntv)
    with torch.no_grad():
        pred_alt = head(rna_hidden, dna_alt)
    assert not torch.equal(pred_alt, pred), "DNA input has no effect on prediction"
    print("PASS  DNA input affects prediction (head is wired through both inputs)")

    # Gradient flow.
    head.train()
    pred = head(rna_hidden, dna_repr)
    pred.pow(2).mean().backward()
    n_grad = sum(1 for p in head.parameters() if p.grad is not None and p.grad.abs().sum() > 0)
    n_total = sum(1 for p in head.parameters())
    assert n_grad == n_total, f"only {n_grad}/{n_total} head params got gradient"
    print(f"PASS  {n_grad}/{n_total} head params received nonzero gradient")

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
