"""Bit-exact day-0 equivalence test for CoreWithRNA.

Verifies that CoreWithRNA is identical to stock NTv3 Core when rna_bias is
None or zeros, and that a nonzero rna_bias actually changes the output (so
we know the line is wired in, not just a no-op).
"""
import os
import sys
from pathlib import Path

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

# Make `sc_nona/model/...` importable when running this file directly.
_REPO = Path(__file__).resolve().parents[2]  # sc_nona/
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.core_with_rna import CoreWithRNA, load_pretrained  # noqa: E402

MODEL = "InstaDeepAI/NTv3_8M_pre"
TOKEN = os.environ.get("HF_TOKEN")


def main() -> None:
    torch.manual_seed(0)

    # Stock NTv3 (downloads upstream modeling code via trust_remote_code).
    stock = AutoModelForMaskedLM.from_pretrained(
        MODEL, trust_remote_code=True, token=TOKEN
    )
    stock.eval()

    # Our subclass loaded from the same weights.
    ours = load_pretrained(MODEL, token=TOKEN)
    ours.eval()

    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True, token=TOKEN)
    seq = "ACGT" * 128  # 512 chars -> 4 bins at L/128
    batch = tok(
        [seq],
        add_special_tokens=False,
        padding=True,
        pad_to_multiple_of=128,
        return_tensors="pt",
    )

    B, L = batch["input_ids"].shape
    bottleneck_L = L // (2 ** ours.config.num_downsamples)
    embed_dim = ours.config.embed_dim
    print(f"B={B}, L={L}, L/128={bottleneck_L}, embed_dim={embed_dim}")

    with torch.no_grad():
        out_stock = stock.core(batch["input_ids"])
        out_none = ours(batch["input_ids"])
        out_zero = ours(
            batch["input_ids"],
            rna_bias=torch.zeros(B, bottleneck_L, embed_dim),
        )
        out_one = ours(
            batch["input_ids"],
            rna_bias=torch.ones(B, bottleneck_L, embed_dim),
        )

    # 1. CoreWithRNA(rna_bias=None) is bit-exact w.r.t. stock NTv3.
    assert torch.equal(out_none["logits"], out_stock["logits"]), (
        "rna_bias=None path diverged from stock NTv3 -- forward body drifted"
    )
    print("PASS  rna_bias=None  is bit-exact-equal to stock NTv3")

    # 2. Zero-init bias path is bit-exact w.r.t. None path.
    assert torch.equal(out_zero["logits"], out_stock["logits"]), (
        "rna_bias=zeros path diverged from stock NTv3"
    )
    print("PASS  rna_bias=zeros is bit-exact-equal to stock NTv3")

    # 3. Nonzero bias must change the output (otherwise the line is dead).
    assert not torch.equal(out_one["logits"], out_stock["logits"]), (
        "rna_bias=ones produced identical output -- the bias is NOT wired in!"
    )
    max_delta = (out_one["logits"] - out_stock["logits"]).abs().max().item()
    print(f"PASS  rna_bias=ones changes output (max|delta logits| = {max_delta:.4f})")

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
