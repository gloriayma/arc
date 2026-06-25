"""Shape + day-0 equivalence + gradient checks for DNARNAHead.

Day-0 equivalence: when initialized via `from_pretrained_lm_head`, the head's
output must equal applying the stock NTv3 lm_head to dna_repr alone (regardless
of rna_stream), because the RNA-input weight slice is zero.
"""
import sys
from pathlib import Path
from types import SimpleNamespace

import torch
from torch import nn
from torch.nn import functional as F

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.dna_rna_head import DNARNAHead  # noqa: E402


def _stock_lm_head_forward(lm_head: nn.ModuleDict, y: torch.Tensor) -> torch.Tensor:
    """Replica of NTv3's Core.lm_head application (lines 822-827 in modeling)."""
    y = F.gelu(y, approximate="tanh")
    for hl in lm_head["hidden_layers"]:
        y = F.gelu(hl(y), approximate="tanh")
    return lm_head["head"](y)


def _build_stock_lm_head(embed_dim, conv_init_embed_dim, num_hidden, alphabet_size):
    return nn.ModuleDict({
        "hidden_layers": nn.ModuleList(
            [nn.Linear(embed_dim, embed_dim) for _ in range(num_hidden)]
        ),
        "head": nn.Linear(conv_init_embed_dim, alphabet_size),
    })


def run_case(num_hidden: int) -> None:
    torch.manual_seed(0)
    B, L_bin = 2, 4
    embed_dim = 256
    conv_init_embed_dim = embed_dim
    alphabet_size = 11
    num_downsamples = 7  # L/128
    d_rna_inject = 64
    L = L_bin * (2 ** num_downsamples)

    cfg = SimpleNamespace(
        embed_dim=embed_dim,
        conv_init_embed_dim=conv_init_embed_dim,
        num_hidden_layers_head=num_hidden,
        alphabet_size=alphabet_size,
        num_downsamples=num_downsamples,
    )

    dna_repr = torch.randn(B, L, embed_dim)
    rna_stream = torch.randn(B, L_bin, d_rna_inject)

    stock = _build_stock_lm_head(embed_dim, conv_init_embed_dim, num_hidden, alphabet_size)
    stock.eval()
    with torch.no_grad():
        stock_logits = _stock_lm_head_forward(stock, dna_repr)

    head = DNARNAHead.from_pretrained_lm_head(cfg, d_rna_inject, stock)
    head.eval()

    # Shape.
    with torch.no_grad():
        logits = head(dna_repr, rna_stream)
    assert logits.shape == (B, L, alphabet_size), (
        f"wrong shape {tuple(logits.shape)} (num_hidden={num_hidden})"
    )

    # Day-0 equivalence: nonzero rna_stream MUST NOT affect output.
    assert torch.allclose(logits, stock_logits, atol=1e-6), (
        f"day-0 logits diverged from stock (num_hidden={num_hidden}); "
        f"max delta = {(logits - stock_logits).abs().max().item()}"
    )

    # Random init (no copy) should NOT match stock.
    head_fresh = DNARNAHead(cfg, d_rna_inject)
    head_fresh.eval()
    with torch.no_grad():
        logits_fresh = head_fresh(dna_repr, rna_stream)
    assert not torch.allclose(logits_fresh, stock_logits, atol=1e-4), (
        "random-init DNARNAHead matches stock -- something's off"
    )

    # After training would move RNA weights off zero, RNA must affect output. Simulate by
    # filling the RNA slice of the first weight with small nonzero values.
    with torch.no_grad():
        if num_hidden > 0:
            head.hidden_layers[0].weight[:, embed_dim:].normal_(std=0.1)
        else:
            head.head.weight[:, conv_init_embed_dim:].normal_(std=0.1)
    with torch.no_grad():
        logits_with_rna = head(dna_repr, rna_stream)
    assert not torch.allclose(logits_with_rna, stock_logits, atol=1e-4), (
        "even with RNA slice nonzero, output unchanged -- RNA is not wired in!"
    )

    # Gradient flow.
    head.train()
    logits = head(dna_repr, rna_stream)
    logits.pow(2).mean().backward()
    n_grad = sum(1 for p in head.parameters() if p.grad is not None and p.grad.abs().sum() > 0)
    n_total = sum(1 for p in head.parameters())
    assert n_grad == n_total, f"only {n_grad}/{n_total} params got gradient"

    print(f"PASS  num_hidden={num_hidden}: shape, day-0 equivalence, RNA wiring, gradient flow")


def main() -> None:
    for num_hidden in (0, 2):
        run_case(num_hidden)
    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
