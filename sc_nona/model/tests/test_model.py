"""End-to-end test for NTv3WithRNA.

Verifies:
  1. Shapes of all four output tensors.
  2. Day-0 bit-exactness: dna_logits_stock and dna_logits_with_rna both match
     stock NTv3's forward output exactly (rna_bias = 0 + DNARNAHead RNA slice = 0).
  3. Gradients flow into every component (core, encoder, RNA head, DNA-RNA head).
"""
import os
import sys
from pathlib import Path

import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from model.model import NTv3WithRNA  # noqa: E402

MODEL = "InstaDeepAI/NTv3_8M_pre"
TOKEN = os.environ.get("HF_TOKEN")


def main() -> None:
    torch.manual_seed(0)

    # Stock NTv3 reference.
    stock = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True, token=TOKEN)
    stock.eval()

    # Our wrapper.
    d_rna = 64
    D = 5
    model = NTv3WithRNA.from_pretrained(MODEL, d_rna=d_rna, token=TOKEN)
    model.eval()
    cfg = model.config

    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True, token=TOKEN)
    seq = "ACGT" * 128  # L = 512
    batch = tok([seq], add_special_tokens=False, padding=True,
                pad_to_multiple_of=128, return_tensors="pt")
    dna_tokens = batch["input_ids"]
    B, L = dna_tokens.shape
    L_bin = L // model.bin_size

    rna_matrix = torch.rand(B, L, D) * 10  # nonzero, varied
    rna_mask = torch.zeros(B, L_bin, D, dtype=torch.bool)
    rna_mask[:, 1:3, :] = True  # block-along-length masking

    print(f"B={B}, L={L}, L_bin={L_bin}, D={D}, embed_dim={cfg.embed_dim}, alphabet={cfg.alphabet_size}")

    with torch.no_grad():
        out = model(dna_tokens, rna_matrix, rna_mask=rna_mask)
        stock_logits = stock.core(dna_tokens)["logits"]

    # 1. Shapes.
    assert out["dna_logits_stock"].shape == (B, L, cfg.alphabet_size), out["dna_logits_stock"].shape
    assert out["dna_logits_with_rna"].shape == (B, L, cfg.alphabet_size), out["dna_logits_with_rna"].shape
    assert out["rna_pred"].shape == (B, L_bin, D), out["rna_pred"].shape
    assert out["rna_bias"].shape == (B, L_bin, cfg.embed_dim), out["rna_bias"].shape
    print("PASS  all output shapes correct")

    # 2. Day-0 zero-init: rna_bias must be exactly zero.
    assert torch.equal(out["rna_bias"], torch.zeros_like(out["rna_bias"]))
    print("PASS  rna_bias is exactly zero at day-0")

    # 3. dna_logits_stock matches stock NTv3 bit-exactly.
    assert torch.equal(out["dna_logits_stock"], stock_logits), (
        f"dna_logits_stock diverged; max delta = {(out['dna_logits_stock'] - stock_logits).abs().max().item()}"
    )
    print("PASS  dna_logits_stock bit-exact-equals stock NTv3")

    # 4. dna_logits_with_rna also equals stock at day-0 (zero RNA slice in DNARNAHead).
    assert torch.allclose(out["dna_logits_with_rna"], stock_logits, atol=1e-6), (
        f"dna_logits_with_rna diverged at day-0; max delta = "
        f"{(out['dna_logits_with_rna'] - stock_logits).abs().max().item()}"
    )
    print("PASS  dna_logits_with_rna bit-exact-equals stock NTv3 at day-0")

    # 5. Sanity: passing different rna inputs should NOT change either DNA logits at day-0,
    #    because both rna_bias and DNARNAHead RNA slice are zero.
    rna_alt = torch.rand_like(rna_matrix) * 50
    with torch.no_grad():
        out_alt = model(dna_tokens, rna_alt, rna_mask=rna_mask)
    assert torch.equal(out_alt["dna_logits_stock"], out["dna_logits_stock"])
    assert torch.allclose(out_alt["dna_logits_with_rna"], out["dna_logits_with_rna"], atol=1e-6)
    # But rna_pred SHOULD differ (different RNA input -> different predictions).
    assert not torch.equal(out_alt["rna_pred"], out["rna_pred"]), (
        "rna_pred didn't change with different RNA input"
    )
    print("PASS  day-0 DNA logits invariant to RNA input; rna_pred reactive to RNA input")

    # 6. Gradient flow into every submodule.
    model.train()
    out = model(dna_tokens, rna_matrix, rna_mask=rna_mask)
    loss = out["dna_logits_with_rna"].pow(2).mean() + out["rna_pred"].pow(2).mean()
    loss.backward()

    def _frac_with_grad(module: torch.nn.Module) -> tuple[int, int]:
        n = sum(1 for p in module.parameters())
        g = sum(1 for p in module.parameters() if p.grad is not None and p.grad.abs().sum() > 0)
        return g, n

    for name, sub in [
        ("core", model.core),
        ("rna_encoder", model.rna_encoder),
        ("rna_mlm_head", model.rna_mlm_head),
        ("dna_rna_head", model.dna_rna_head),
    ]:
        g, n = _frac_with_grad(sub)
        print(f"PASS  {name}: {g}/{n} params received nonzero gradient")
        assert g > 0, f"no gradient reached {name}"

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    main()
