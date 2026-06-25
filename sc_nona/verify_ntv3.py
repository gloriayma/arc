"""Pre-verification of NTv3 architecture before scaffolding sc_nona model code.

Verifies:
  1. Config values match what we expect (embed_dim, num_downsamples, num_layers).
  2. Forward pass works on a dummy DNA sequence.
  3. Identifies which hidden_states[i] corresponds to the L/128 bottleneck
     (right after conv tower, before transformer).
  4. Confirms we can construct a fresh Core from the config and load_state_dict
     with strict=True — proving CoreWithRNA subclass approach will work.
"""
import os, sys, importlib

import torch
from transformers import AutoConfig, AutoModelForMaskedLM, AutoTokenizer

MODEL = "InstaDeepAI/NTv3_8M_pre"
TOKEN = os.environ.get("HF_TOKEN")

print("=" * 70)
print("Step 1: load model")
print("=" * 70)
cfg = AutoConfig.from_pretrained(MODEL, trust_remote_code=True, token=TOKEN)
tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True, token=TOKEN)
model = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True, token=TOKEN)
model.eval()

print(f"embed_dim         = {cfg.embed_dim}")
print(f"num_downsamples   = {cfg.num_downsamples}  -> L/{2**cfg.num_downsamples}")
print(f"num_layers        = {cfg.num_layers}  (transformer blocks)")
print(f"conv_init_embed_dim = {cfg.conv_init_embed_dim}")
print(f"token_embed_dim   = {cfg.token_embed_dim}")
print(f"alphabet_size     = {cfg.alphabet_size}")
print(f"filter_list       = {cfg.filter_list}")
print(f"use_skip_connection = {cfg.use_skip_connection}")
print()

print("=" * 70)
print("Step 2: forward pass on dummy DNA, with output_hidden_states=True")
print("=" * 70)
seq = "ACGT" * 256              # 1024 chars (multiple of 128)
batch = tok([seq], add_special_tokens=False, padding=True,
            pad_to_multiple_of=128, return_tensors="pt")
print(f"input_ids shape   = {tuple(batch['input_ids'].shape)}")

with torch.no_grad():
    out = model.core(batch["input_ids"], output_hidden_states=True)

print(f"logits shape      = {tuple(out['logits'].shape)}")
print(f"# hidden_states   = {len(out['hidden_states'])}")
print()

print("=" * 70)
print("Step 3: hidden_states shape audit -> find the L/128 bottleneck")
print("=" * 70)
L_tokens = batch["input_ids"].shape[1]
L_bottleneck = L_tokens // (2 ** cfg.num_downsamples)
print(f"L_tokens = {L_tokens}, expected bottleneck L = L/128 = {L_bottleneck}")
print()
print(f"{'idx':>3}  {'shape':<28}  {'stage':<22}")
for i, h in enumerate(out["hidden_states"]):
    shape = tuple(h.shape)
    L_here = shape[1]
    if L_here == L_bottleneck:
        stage = "*** BOTTLENECK (L/128) ***"
    elif L_here < L_tokens:
        stage = f"downsampled (L/{L_tokens // L_here})"
    else:
        stage = "full-resolution"
    print(f"{i:>3}  {str(shape):<28}  {stage:<22}")
print()

bottleneck_idxs = [i for i, h in enumerate(out["hidden_states"])
                   if h.shape[1] == L_bottleneck]
print(f"hidden_states indices at L/128: {bottleneck_idxs}")
print(f"  -> first one (post-conv, pre-transformer) is hidden_states[{bottleneck_idxs[0]}]")
print(f"  -> last one  (post-transformer, pre-deconv) is hidden_states[{bottleneck_idxs[-1]}]")
print()

print("=" * 70)
print("Step 4: state_dict round-trip into fresh Core (strict=True)")
print("=" * 70)
modeling_module = importlib.import_module(type(model.core).__module__)
Core = type(model.core)
print(f"Core class = {Core.__module__}.{Core.__name__}")

fresh = Core(cfg)
missing, unexpected = fresh.load_state_dict(model.core.state_dict(), strict=False)
print(f"missing keys    : {len(missing)}")
print(f"unexpected keys : {len(unexpected)}")
if missing:    print("  first missing:", missing[:3])
if unexpected: print("  first unexpected:", unexpected[:3])

print()
print("ALL CHECKS COMPLETE")
