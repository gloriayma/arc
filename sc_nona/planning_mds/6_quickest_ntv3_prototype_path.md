# Quickest Path to a Working NTv3 + scRNA Prototype

**Audience:** Gloria Ma (gloria.ma@arcinstitute.org), Arc Institute
**Date:** 2026-06-17
**Scope:** A concrete, hands-on plan to fine-tune **pretrained** NTv3 (Nucleotide Transformer v3) with scRNA-seq coverage as a conditioning signal, in a way that is as drop-in as possible with the HuggingFace `transformers` API. This document supersedes nothing — it is a tactical companion to `5_architectures_for_rna_conditioned_dna_lms.md`, which surveys the architectural space. The recommendations here are calibrated against what NTv3 actually exposes on HuggingFace as of late February 2026.

---

## TL;DR

1. **Use `InstaDeepAI/NTv3_100M_pre` as the starting checkpoint.** 100M is the sweet spot for fast iteration. Drop down to `NTv3_8M_pre` for smoke tests, scale up to `NTv3_650M_pre` once the recipe works. All loaded with `AutoModelForMaskedLM.from_pretrained(..., trust_remote_code=True)`.
2. **The single most drop-in conditioning method for NTv3 is per-token additive embedding on the input.** NTv3 is single-base char-level (vocab=11, embed dim=16, then projected to 768/1536). The bottleneck of injecting one extra signal into the embedding stream is ~10 lines of code. This is approach (b) in your question — it works because NTv3 already does its first big projection inside the conv stem; a small additive bias before the conv stem is benign.
3. **The second-most drop-in pattern, and probably better long-term, is to reuse NTv3's existing FiLM-style adaptive-layer-norm conditioning channel.** NTv3 has *baked-in* adaptive layer norm machinery for species conditioning (see `adaptive_layers.py` in the InstaDeep repo). Repurposing this for a learned RNA condition is architecturally identical to what NTv3 was already trained to handle — no surgery beyond providing a different condition vector.
4. **Skip the "posttrained" checkpoint.** Posttrained NTv3 is supervised on ~16k functional tracks (BigWig/BED). Starting from it means you fight the existing track heads and you bake in a head shape that doesn't match scRNA Smart-seq coverage. Starting from `_pre` is correct.
5. **Compute target for prototype: one H100 or A100 80 GB, 100M params, ~32 kb windows.** Trainable in days, not weeks.
6. **Recommended training objective for the smoke test: masked-nucleotide-LM with RNA-coverage conditioning + a separate auxiliary coverage-reconstruction head.** The MLM loss validates the model still does its job; the auxiliary head validates the conditioning channel actually carries information.
7. **The fastest, most credible alternative path is scooby (Borzoi+LoRA at single-cell resolution).** It exists, it works, and it solves a strict subset of your problem (RNA as output, not input). Calibrate against scooby — it is the strongest baseline.

---

## 1. NTv3 Checkpoints on HuggingFace — Exact Strings, Sizes, Specs

NTv3 dropped on 2025-12-23 (paper: bioRxiv 2025.12.22.695963; press release at instadeep.com/2026/02/modelling-the-genome-with-ntv3/). The full HuggingFace collection lives at https://huggingface.co/collections/InstaDeepAI/nucleotide-transformer-v3.

### 1.1 Pretrained ("_pre"): the right starting point

Pretrained checkpoints are trained with **masked LM over 9T base pairs from OpenGenome2**, covering >128,000 species. They expose only the MLM logits + hidden states. Use these.

| HF Model ID | Params | Hidden dim | FFN dim | Heads | Transformer layers | Downsample stages | Notes |
|---|---|---|---|---|---|---|---|
| `InstaDeepAI/NTv3_8M_pre` | 7.69M | 256 | 1024 | 8 | 2 | 7 | Smallest. Use for smoke tests + unit tests. |
| `InstaDeepAI/NTv3_100M_pre` | 100M | 768 | 3072 | 12 | 6 | 7 | **Recommended starting point.** Best size/performance trade-off. |
| `InstaDeepAI/NTv3_650M_pre` | 700M | 1536 | 6144 | 24 | 12 | 7 | Largest. Use after the 100M recipe works. |
| `InstaDeepAI/NTv3_5downsample_pre` | ~600M | ? | ? | ? | ? | 5 | 5 downsamples → 32× compression instead of 128×. Requires divisible-by-32 input. Higher resolution in the transformer interior. |

**8kb-only variants** (do NOT use as base — these are length-restricted exploration models):
- `InstaDeepAI/NTv3_8M_pre_8kb`
- `InstaDeepAI/NTv3_100M_pre_8kb`
- `InstaDeepAI/NTv3_650M_pre_8kb`
- `InstaDeepAI/NTv3_5downsample_pre_8kb`

These are the curriculum stage-1 checkpoints (trained on 1–8 kb windows). The full-context checkpoints went through stage-1 plus length-extension to 1 Mb. The base "_pre" models are the right thing for a 32–131 kb prototype window.

### 1.2 Posttrained ("_post") — skip these

Posttrained checkpoints add task-specific heads for ~16k BigWig/BED tracks across 24 species and adaptive species conditioning. They return `bigwig_tracks_logits` and `bed_tracks_logits` in addition to `logits`. Quoting the model card for `NTv3_100M_post`:

> "This is a **post-trained (supervised) multi-species model** that can predict functional genomics tracks (BigWig) and genomic elements (BED) across multiple species. It builds on the pre-trained NTv3 model with additional conditioning mechanisms and task-specific heads."

Two reasons to skip posttrained for your project:
1. You're injecting RNA as **input**, not predicting RNA as a separate track. The posttrained track heads are a different output shape than scRNA coverage and you'd just be deleting them.
2. The posttrained model is conditioned on *species* via adaptive-layer-norm. That conditioning channel is exactly the one you want to commandeer for RNA. Starting from posttrained means fighting the species conditioning; starting from pre means writing on a clean slate (the adaptive layer norms are present in code but the parameters are uninitialized / pass-through in the pre checkpoint).

| HF Model ID | Params | Notes |
|---|---|---|
| `InstaDeepAI/NTv3_100M_post` | 100M | Posttrained, species-conditioned, has track heads. |
| `InstaDeepAI/NTv3_650M_post` | 700M | Same as above, 650M scale. |
| `InstaDeepAI/NTv3_100M_post_131kb` | 100M | 131 kb context variant. |
| `InstaDeepAI/NTv3_650M_post_131kb` | 700M | 131 kb context variant. |
| `InstaDeepAI/NTv3_5downsample_post` | ~600M | 5-downsample posttrained. |
| `InstaDeepAI/NTv3_5downsample_post_131kb` | ~600M | 5-downsample, 131 kb. |

### 1.3 Other variants

- `InstaDeepAI/NTv3_generative` (0.7B): a posttrained model further fine-tuned with **masked-diffusion language modeling (MDLM)**, intended for generating sequences with target activity. Out of scope.
- `InstaDeepAI/ntv3_base_model`: a base model checkpoint, likely an alias.

### 1.4 Quick load test (copy-pasteable)

```python
from transformers import AutoTokenizer, AutoModelForMaskedLM

repo = "InstaDeepAI/NTv3_100M_pre"
tok   = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained(repo, trust_remote_code=True)

seq = "ATCG" * 32  # 128 bp; must be multiple of 128 = 2**num_downsamples
batch = tok([seq], add_special_tokens=False,
            padding=True, pad_to_multiple_of=128, return_tensors="pt")

out = model(**batch, output_hidden_states=True)
print(out.logits.shape)            # (1, 128, 11)
print(out.hidden_states[-1].shape) # (1, 128, 768) for 100M
```

> **Gating note:** the InstaDeepAI HF repos require accepting the license. Configure `HF_HUB_TOKEN` and visit each model page once to accept terms. The custom modeling file (e.g. `ntv3_huggingface_new.py`) is bundled inside the repo and loaded via `trust_remote_code=True`.

---

## 2. HuggingFace Interface — Exactly How NTv3 Plugs In

### 2.1 The classes you need

- `AutoTokenizer` (character-level over `A T C G N` + `<unk> <pad> <mask> <cls> <eos> <bos>`).
- `AutoModelForMaskedLM` — primary entry. Returns `MaskedLMOutput`-like object with `logits`, `hidden_states`, `attentions`, and (for posttrained) extra fields.
- `AutoModel` — same network but without the LM head. Returns embeddings.

Both rely on `trust_remote_code=True`. The custom config class is `Ntv3PreTrainedConfig` and is importable as:

```python
from ntv3_huggingface_new import Ntv3PreTrainedConfig
```

### 2.2 Forward signature

For the **pretrained** variant:

```python
out = model(
    input_ids=batch["input_ids"],            # (B, L), L % 128 == 0
    attention_mask=batch["attention_mask"],  # standard HF
    output_hidden_states=True,               # for adapter insertion / probing
    output_attentions=True,                  # for interpretability
    labels=None,                             # MLM labels if training
)

# Outputs:
out.logits           # (B, L, 11)
out.hidden_states    # tuple of tensors per layer, all at full L resolution
out.attentions       # tuple of attn maps from each transformer layer (at compressed L/128)
out.loss             # cross-entropy when labels passed
```

For the **posttrained** variant the same signature plus:

```python
out = model(
    input_ids=...,
    species_ids=species_ids,   # required
)
out.bigwig_tracks_logits      # (B, L_cropped, num_tracks)
out.bed_tracks_logits         # (B, L_cropped, num_classes)
out.embedding                 # post-deconv tower (B, L, D)
out.after_transformer_embedding   # post-transformer (compressed, B, L/128, D)
```

`L_cropped` is the middle 37.5% of L (Borzoi-style crop).

### 2.3 Where you can hook in (per the model card, quoted)

The 100M_pre card explicitly documents these hook points:

```python
config = Ntv3PreTrainedConfig.from_pretrained(repo)
config.embeddings_layers_to_save = (1, 2)             # save outputs of these transformer layers
config.attention_maps_to_save = [(1, 0), (2, 1)]      # (layer, head) pairs
config.deconv_layers_to_save = (1, 2)                 # save deconv-tower intermediate outputs

model = AutoModelForMaskedLM.from_pretrained(repo, config=config, trust_remote_code=True)
core_out = model.core(**batch, output_hidden_states=True, output_attentions=True)
```

The fact that there is a `model.core` attribute exposing intermediate-state APIs is good news. It means hooks (`register_forward_hook`) work cleanly and you can inject without monkey-patching the inner modules.

### 2.4 The architecture, in code-order

```
input_ids (B, L)             # L must be multiple of 128
    |
    v
Embedding (V=11 -> embed_dim=16)
    |
    v
Stem: 1-D conv, kernel=15, projects to model_dim (256/768/1536)
    |
    v
ConvTower (7 downsamples):
    [ConvBlock -> ResidualConvBlock -> avg_pool(2)] x 7
    Each downsample halves L; total reduction = 128x.
    Residuals saved at each level for skip connections.
    |
    v
Transformer stack (2 / 6 / 12 layers depending on size):
    SelfAttentionBlock with:
      - Multi-head self-attention
      - Rotary positional embeddings (RoPE, UPPER_FREQ=10000, optional Yarn rescaling)
      - Pre-layer-norm
      - GLU-style FFN (optional)
    |
    v
DeconvTower (7 upsamples, mirrors ConvTower):
    [DeConvBlock -> ResidualDeConvBlock] with skip connection: x = x + residuals[level]
    |
    v
LMHead: stacked Linear+GELU layers -> (B, L, V=11)
```

Adaptive-layer-norm variants (`AdaptiveLayerNorm`, `AdaptiveSelfAttentionBlock`, `ConditionedConvTowerBlock`, `ConditionedDeConvTowerBlock`) are defined in `adaptive_layers.py` and used in the posttrained variant for species conditioning. They implement FiLM: `output = x * (1 + gamma(c)) + beta(c)` where `c` is a condition embedding. **These layers are present in the codebase but inert (or absent) in `_pre` checkpoints** — they're the natural injection point.

### 2.5 What's NOT specified in the public model card

- The exact PyTorch class names of the conv blocks (the HF code is auto-generated from a JAX base via Equinox/nnx → PyTorch).
- Whether layer norms freeze cleanly under LoRA without numerical surprises. Expect this is fine but verify with a tiny smoke test.
- Whether `attention_mask` is properly propagated through the conv + transformer stack. NTv3's attention operates on the *compressed* sequence; the mask is downsampled implicitly. **This matters for variable-length batches.** Easiest: pad to a fixed window length per batch (32k or 65k) and forget about attention_mask edge cases.

### 2.6 Official notebooks to crib from

In `github.com/instadeepai/nucleotide-transformer/tree/main/notebooks/nucleotide_transformer_v3`:
- `inference_pretrained.ipynb` — loading + forward pass for `_pre`.
- `inference_posttrained.ipynb` — loading + species conditioning + track readout for `_post`.

There is **no official fine-tuning notebook for NTv3 yet** (as of February 2026 — checked Feb 2026 dates on the HF collection). Older NTv2 fine-tuning notebooks exist (LoRA, PEFT) and the patterns transfer.

---

## 3. NTv3 Context Length & Token Resolution — Critical for RNA Alignment

### 3.1 Effective context

| Variant | Pretrained input window | Effective max | Notes |
|---|---|---|---|
| `_pre` (default) | 1 Mb | 1 Mb at base resolution | Curriculum: stage-1 8 kb → stage-2 length-extension to 1 Mb. |
| `_pre_8kb` | 8 kb | 8 kb | Stage-1 only. |
| `_post_131kb` | 131 kb | 131 kb | Posttrained at this length explicitly. |

For a prototype, **start with 32 kb windows.** The NTv3 paper's downstream benchmark is at 32 kb; the model is well-behaved at that length; 32 kb covers a typical gene + nearby regulatory landscape; memory is tractable.

If you want enhancer-promoter pairings at TAD scale, **bump to 131 kb** (use `_pre` model, NOT the 131kb-only variant) once the 32 kb pipeline works.

### 3.2 Token resolution

- **Input tokenization:** 1 bp per token. Vocab=11.
- **Transformer interior:** 1 token = 128 bp after 7 downsamples (or 32 bp for the 5-downsample variant).
- **Output (after deconv):** back to 1 bp per token via skip connections.

So you have two natural resolutions to align RNA coverage to:

| Resolution | Where in the model | Tokens per 32 kb window | Tokens per 131 kb | Tokens per 1 Mb |
|---|---|---|---|---|
| 1 bp | input + final output | 32,768 | 131,072 | 1,048,576 |
| 128 bp | transformer interior | 256 | 1,024 | 8,192 |
| 32 bp | 5-downsample interior | 1,024 | 4,096 | 32,768 |

**This is a beautiful alignment with Enformer/Borzoi conventions.** Enformer uses 128 bp bins. Borzoi predicts at 32 bp bins. So:

- **If you bin coverage to 128 bp**, inject at the transformer interior (post-conv-tower, pre-transformer). The number of RNA tokens equals the number of DNA tokens. Cleanest math.
- **If you bin coverage to 1 bp**, inject at input or pre-stem. Largest dimensionality but cleanest semantics.
- **If you bin to 32 bp**, use the 5-downsample variant.

**My recommendation: bin coverage to 128 bp and inject at the transformer interior.** This matches the model's natural rate and matches the conventions of all the supervised expression models you'd compare against.

### 3.3 Coverage-binning scheme (what aggregation function?)

The user's CLAUDE.md asks specifically about BigWig-style coverage at full bp resolution. To get from per-bp coverage to 128 bp tokens:

| Aggregation | Pros | Cons | Recommended? |
|---|---|---|---|
| **Mean** of bp coverage in bin | Standard for Enformer/Borzoi. Interpretable. | Dilutes sharp signal (splice sites, TSS spikes). | Default. |
| **Sum** | Preserves total read mass. | Larger dynamic range. | Use if you log-transform downstream anyway. |
| **Max** | Preserves sharp peaks (TSS, splice sites). | Drops contextual mean. | Useful as a second channel. |
| **Log1p(mean)** | Compresses the dynamic range; standard for scRNA. | Loses sign of any zero/nonzero distinction at low coverage. | **Recommended primary signal.** |
| **Multi-channel** (log1p-mean + log1p-max) | Captures both. | 2× dimensionality. | Recommended for v0.2. |

For Smart-seq full-length scRNA, log1p(mean coverage) at 128 bp resolution per cell is the canonical signal. Stack across cells → `(num_cells, L/128)` matrix per locus.

---

## 4. Minimum-Surgery Integration Patterns — Ranked

This is the meat of your question. I'll rank from most-drop-in to most-invasive, with concrete code sketches and effort estimates.

### Approach A — Soft-prompt prefix tokens (the lightest possible touch)

**Idea.** Compress the RNA tracks for this locus into K learned "context tokens" (say K=64). Prepend to the DNA sequence. NTv3 attends to them through standard self-attention. Zero changes to NTv3's internals.

**Problem with NTv3 specifically.** NTv3's first 7 layers are convolutions, not attention. A prefix token concatenated with the DNA input will be convolved with the rest of the sequence and lose its meaning as a "prefix." This pattern is well-suited to a pure-transformer LM, not a conv-front-end LM. **Not recommended for NTv3.**

You *could* concatenate prefix tokens after the conv tower (at the transformer interior), and it works mathematically. But at that point you've already done architectural surgery, and approach C below is cleaner.

**Effort:** 30 LOC if it worked. ~0 LOC saved over approach C in practice.

**Verdict:** Skip.

---

### Approach B — Per-token additive embedding on the input (the "Enformer-style track injection")

**Idea.** For each input nucleotide token, you have RNA coverage values at that position. Embed those values (a learned MLP from a fixed-dim coverage vector to NTv3's input embed dim, which is 16), and add to the token embedding before NTv3 sees it.

**Why this is genuinely drop-in.** NTv3's input pipeline is:
```
input_ids -> Embedding(V=11, D=16) -> Stem(conv k=15, D=16 -> D=256/768/1536)
```
We add our RNA embedding right after the `Embedding` and before the `Stem`. From the rest of the network's perspective, the input embedding distribution is now shifted by a learned, RNA-conditioned offset. The Stem conv smooths over it. Zero changes inside NTv3.

**Concrete code.** Wrap NTv3 with a thin `nn.Module`:

```python
import torch
import torch.nn as nn
from transformers import AutoModelForMaskedLM

class NTv3WithRNAInput(nn.Module):
    def __init__(self, repo="InstaDeepAI/NTv3_100M_pre", rna_channels=1):
        super().__init__()
        self.ntv3 = AutoModelForMaskedLM.from_pretrained(repo, trust_remote_code=True)
        # NTv3 token embed dim = 16 (small, before stem-conv expands to 256/768/1536)
        embed_dim = self.ntv3.get_input_embeddings().embedding_dim
        self.rna_proj = nn.Sequential(
            nn.Linear(rna_channels, 32),
            nn.GELU(),
            nn.Linear(32, embed_dim),
        )
        # Zero-init the final layer so initial behavior is identical to vanilla NTv3
        nn.init.zeros_(self.rna_proj[-1].weight)
        nn.init.zeros_(self.rna_proj[-1].bias)

    def forward(self, input_ids, rna_coverage, attention_mask=None, labels=None,
                output_hidden_states=False):
        # rna_coverage: (B, L, C) where C = number of RNA channels (e.g., 1 = log1p mean cov)
        # Compute the token embeddings explicitly, add RNA, then pass via inputs_embeds.
        tok_emb = self.ntv3.get_input_embeddings()(input_ids)   # (B, L, 16)
        rna_emb = self.rna_proj(rna_coverage)                    # (B, L, 16)
        inputs_embeds = tok_emb + rna_emb
        return self.ntv3(
            inputs_embeds=inputs_embeds,
            attention_mask=attention_mask,
            labels=labels,
            output_hidden_states=output_hidden_states,
        )
```

> **Caveat.** The above assumes `inputs_embeds` is a supported keyword in NTv3's HF wrapper. Many custom-code HF models implement it; some don't. **First verification step**: read `ntv3_huggingface_new.py` and check the `forward()` signature. If `inputs_embeds` isn't supported, you have two fallbacks:
> 1. Use `register_forward_pre_hook` on the input `nn.Embedding` to add the RNA embedding.
> 2. Replace `model.embedding` with a wrapped module that adds RNA before returning.
>
> Hooks are easier; I show the hook variant below.

**Hook-based variant (no assumption about `inputs_embeds`):**

```python
import torch
import torch.nn as nn
from transformers import AutoModelForMaskedLM

class RNAConditioner(nn.Module):
    """Stores the current minibatch's RNA tensor and the projection to embed_dim."""
    def __init__(self, rna_channels, embed_dim):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(rna_channels, 32), nn.GELU(), nn.Linear(32, embed_dim),
        )
        nn.init.zeros_(self.proj[-1].weight); nn.init.zeros_(self.proj[-1].bias)
        self.current_rna = None  # set by the dataloader/training loop

    def set(self, rna): self.current_rna = rna
    def __call__(self):  return self.proj(self.current_rna)

def attach_rna_conditioner(ntv3_model, rna_channels=1):
    embed_dim = ntv3_model.get_input_embeddings().embedding_dim
    cond = RNAConditioner(rna_channels, embed_dim).to(next(ntv3_model.parameters()).device)
    def post_embed_hook(module, inputs, output):
        # output: (B, L, embed_dim); add RNA
        if cond.current_rna is None:
            return output
        return output + cond()
    handle = ntv3_model.get_input_embeddings().register_forward_hook(post_embed_hook)
    return cond, handle

# Usage
ntv3 = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/NTv3_100M_pre", trust_remote_code=True)
cond, handle = attach_rna_conditioner(ntv3, rna_channels=1)
# In training loop:
cond.set(rna_tensor_for_this_batch)  # (B, L, 1)
out = ntv3(input_ids=ids)
```

**Pros.** Genuinely the smallest possible code footprint. Survives upstream changes to the model file because you don't touch internals. Zero-init ensures the conditioned model starts as a no-op equivalent of vanilla NTv3 — your loss should not regress on day 1.

**Cons.** Conditioning at the embedding only is a *weak* signal — it has to propagate through 7 conv layers + N transformer layers + 7 deconv layers. The model has many opportunities to ignore it. Mitigations: combine with approach C (cross-attention) or D (full fine-tune).

**Effort.** ~50–100 LOC including the dataloader for RNA, a training loop, and a forward-pass smoke test. **1–2 days of focused work.**

**Risk.** Low.

**Verdict: this is the right starting point. Get this working first.**

---

### Approach C — Repurpose NTv3's existing FiLM / adaptive-layer-norm channel

**Idea.** NTv3's source code defines `AdaptiveLayerNorm`, `AdaptiveSelfAttentionBlock`, `ConditionedConvTowerBlock`, etc. The posttrained models use these to FiLM-condition the network on species. You can either:
1. Replace the species conditioning with RNA conditioning (more invasive — load the post checkpoint, swap in your own condition module).
2. Add your own adaptive-layer-norm modules into the pretrained model.

The cleanest variant of (2): wrap each existing `LayerNorm` in the transformer interior with an additive (FiLM) shift driven by a per-token RNA condition vector. The shift starts at zero (so initial behavior is unchanged) and is learned end-to-end.

**Why this is structurally aligned with NTv3.** The model's authors *already designed it to accept FiLM conditioning*. The "post" checkpoint commits the adaptive-LN parameters to species conditioning. Starting from "pre", those parameters either don't exist (in which case you add them fresh) or are inert. Either way, the inductive bias of FiLM-style conditioning is what the model was trained against during posttraining — so you're working with the grain.

**Concrete approach.** Use PyTorch hooks to wrap each transformer-block layer norm in:

```python
class FiLMWrap(nn.Module):
    """Wrap a LayerNorm to add a learned shift conditioned on RNA at this resolution."""
    def __init__(self, base_ln, model_dim, rna_channels):
        super().__init__()
        self.base_ln = base_ln
        self.model_dim = model_dim
        # Maps per-token RNA to (scale, shift) of model_dim
        self.cond_proj = nn.Linear(rna_channels, 2 * model_dim)
        nn.init.zeros_(self.cond_proj.weight); nn.init.zeros_(self.cond_proj.bias)
        self.rna_at_this_resolution = None  # set externally

    def forward(self, x):
        # x: (B, L_compressed, model_dim)
        out = self.base_ln(x)
        if self.rna_at_this_resolution is None:
            return out
        # rna_at_this_resolution should be (B, L_compressed, rna_channels) — pre-binned
        gs = self.cond_proj(self.rna_at_this_resolution)
        gamma, beta = gs.chunk(2, dim=-1)
        return out * (1 + gamma) + beta
```

Inject by walking the model and replacing every `LayerNorm` in the transformer with `FiLMWrap`. The transformer operates at 128 bp resolution, so you bin RNA coverage at 128 bp once per batch.

**Effort.** ~150–250 LOC. **2–4 days.** The trickiest part is the model-walking code: you have to find the right `LayerNorm` instances and identify which transformer layer each belongs to so you can route the right (B, L_compressed, C) tensor.

**Pros.** Conditions at every block (not just input). Much stronger signal than approach B. Architecturally consistent with NTv3's design.

**Cons.** More code. Requires you to know the model's module structure. Slightly more brittle.

**Verdict.** Excellent v0.2 after approach B works. Don't do this first.

---

### Approach D — Flamingo-style cross-attention adapters (the survey doc's Tier 2 recommendation)

**Idea.** Build a separate RNA encoder (small Hyena/Mamba/Transformer over the (M_cells, L_bins) matrix). Get RNA latent tokens. Insert new cross-attention layers into NTv3 after every Kth transformer block. Freeze NTv3, train only the RNA encoder + cross-attn adapters + LoRA in NTv3.

**This is the *correct* architecture per the survey doc.** It does the job of carrying genuine RNA information into every depth of the DNA model, and it's the canonical "condition a frozen LM on a new modality" pattern.

**Why it's not the v0 prototype.** It requires:
- A bespoke RNA encoder (~300–500 LOC, plus a chunk of experimentation).
- New cross-attention modules (~150 LOC each, inserted into the model).
- A way to interleave the new modules with NTv3's existing forward pass — either monkey-patching or fork.
- Substantially more code review for correctness.

You can absolutely build this — but only after approach B has validated that there's signal in the data and the pipeline works end-to-end.

**Effort.** ~1500–2500 LOC. **2–3 weeks** to a working version that survives ablations.

**Risk.** Medium. The risk is not "it doesn't work" — it's "you spend two weeks before discovering your RNA data is too sparse to be useful."

**Verdict.** v1.0. Build only after approach B shows positive signal.

---

### Approach E — Full continued pretraining

**Idea.** Don't freeze anything; add RNA as a parallel token stream (or as an aux input via approach B) and just continue MLM-style pretraining on the joint (DNA, RNA) dataset.

**Pros.** Maximum flexibility.

**Cons.** Catastrophic forgetting of NTv3's pretraining. Expensive. You're now in 1000s of H100-hours territory, not single-GPU territory.

**Verdict.** Only if approach C/D shows good signal but a frozen backbone is a hard ceiling.

---

### Approach ranking

| Approach | Effort | Risk | When to use |
|---|---|---|---|
| **B. Additive input embedding (zero-init)** | 1–2 days | Low | **v0.1 / smoke test.** Start here. |
| **C. FiLM via adaptive layer norms** | 2–4 days | Low-Med | **v0.2.** Once B shows the pipeline works. |
| **D. Cross-attention adapters** | 2–3 weeks | Med | **v1.0.** Once you've justified deeper conditioning. |
| A. Soft prompt prefix | N/A | High | **Skip** — bad fit for NTv3's conv-front-end. |
| E. Full continued pretraining | Many weeks | High | **Skip for prototype**, revisit after v1.0. |

---

## 5. "Pretrained" vs "Posttrained" — What You're Skipping

Quoting the bioRxiv paper abstract (via the InstaDeep summary):

> "NTv3 is pretrained on 9 trillion base pairs from OpenGenome2 using base-resolution masked language modeling, followed by post-training with a joint objective that integrates continued self-supervision with supervised learning on ∼16,000 functional tracks and annotation labels from 24 animal and plant species."

So:
- **Pre** = pure MLM, no supervised signal. Sequence-only.
- **Post** = same backbone + ~16k functional tracks (BigWig + BED) supervised + species-conditional adaptive layer norms. Returns track logits and annotation logits.

What does starting from "pre" cost you?
1. You skip the supervised functional knowledge (CAGE, ATAC, ChIP, RNA-seq tracks across 24 species). Some of this knowledge would help your task. But:
2. The supervised signal is at the bulk/cell-line level — not single-cell. It's a *different distribution* of expression data than Smart-seq scRNA.
3. The posttrained model's prediction heads are bolted-on and irrelevant to your input-conditioning task.

What does starting from "pre" buy you?
1. A clean slate for your own task heads.
2. No species-conditioning baggage. (Crucial — you want to *override* that channel.)
3. The model is purely an unconditional sequence model. Your conditioning is the first conditioning it sees.

**Verdict: pretrained is the right call.** The user's intuition is correct.

> Caveat: you might benefit from a hybrid — initialize with the post checkpoint, *delete* the track heads, and keep the species adaptive-LN parameters but re-purpose them for cell-type or RNA conditioning. This is approach C taken further. It's a v0.3 experiment, not a v0.1 starting point.

---

## 6. Aligning 128 bp Coverage Bins with NTv3 Tokens

(Mostly covered in §3 — restating the recommendation here.)

- **NTv3 input resolution:** 1 bp/token.
- **NTv3 transformer-interior resolution:** 128 bp/token (after 7 downsamples).
- **Enformer convention:** 128 bp bins.
- **Borzoi convention:** 32 bp bins.

For your prototype:
- **Bin Smart-seq scRNA coverage to 128 bp using log1p(mean).**
- **Inject at the transformer interior** if using approach C/D. The number of RNA tokens (L/128) exactly matches the number of DNA tokens at that level.
- **Inject at input level** (1 bp) if using approach B — but then you must *upsample* the 128 bp bins back to 1 bp (e.g., `repeat_interleave(128, dim=1)`) or *interpolate*. Repeat-interleave is simpler and the conv stem will smooth it.

A double-resolution scheme is also reasonable:
- 1 bp input embedding shift from a high-resolution single-bp coverage (mean over a small kernel).
- 128 bp FiLM at the transformer interior from coarse coverage.

---

## 7. Training Objective for the Prototype

The user listed three candidates. My analysis:

### 7.1 Option A — MLM with RNA conditioning

**The smoke test that proves the wiring works.** Compute MLM cross-entropy with and without the RNA conditioning channel active. If conditioning improves MLM loss on **held-out windows**, your conditioning is carrying information.

**Watch out:** RNA coverage is *correlated* with the underlying DNA (CG content, transcription factor motifs, etc.). Trivial leakage is possible — e.g., the model can use RNA coverage to identify exon vs intron and thus guess masked nucleotides better. That's not a bug; it's actually a valid signal! But you should compare against:
1. **Shuffled-RNA baseline:** same RNA values, randomly permuted across cells. Should *not* improve MLM if your model is using RNA properly.
2. **Reference-aware baseline:** RNA values for a *different* genotype at the same locus. If your model truly learns sequence-RNA covariation, conditional MLM with the right genotype's RNA should outperform conditional MLM with another genotype's RNA.

**Effort:** ~1 day after approach B is implemented.

**Pros:** Cheapest possible validation that conditioning works.

**Cons:** It validates *information flow* but not *biological signal*. A model could improve MLM by latching onto trivial coverage-CG correlation.

### 7.2 Option B — Genotype-conditioned expression prediction

**The serious objective.** Mask one cell's RNA track in a window; predict it from (DNA sequence at this individual's genotype) + (other cells' RNA tracks at this window). Loss is e.g. Poisson NLL or MSE on log1p coverage.

This is the *direct* version of what you want to learn: covariation across cells at the same locus.

**Effort:** ~3–5 days. Requires:
- A masked-cell collation that picks one cell as target, others as input.
- A small prediction head (linear from final hidden state → coverage prediction at 128 bp resolution).
- A Poisson or zero-inflated negative binomial loss.

**Pros:** This is *the* meaningful objective.

**Cons:** More moving parts. Harder to debug if it doesn't work — is the issue in the model, the data, or the loss?

### 7.3 Option C — Promoter-enhancer pairing benchmark

Defer. It's an evaluation, not a training objective. Run it after the model is trained.

### 7.4 Recommendation

**Stage 0:** Approach B (additive input embedding) with **MLM loss only**, no RNA conditioning. This verifies you can load NTv3, run a forward pass, and continue training without crashing. ~1 day.

**Stage 1:** Approach B + **MLM with RNA conditioning** vs **MLM without**. Compare on held-out windows. Smoke-test for information flow. ~1 day.

**Stage 2:** Approach B + **joint MLM loss + auxiliary coverage-reconstruction loss** (Option B). Auxiliary head: linear from `out.hidden_states[-1]` (B, L, D) to (B, L/128, 1), trained with Poisson NLL on the target cell's coverage. The MLM loss keeps NTv3 from degrading; the aux loss validates the conditioning channel. ~3 days.

**Stage 3:** Add genotype-aware evaluation: predict RNA from sequence + leave-one-out cells, evaluate the held-out genotype. ~2 days.

---

## 8. Compute Footprint — Real Numbers

NTv3 isn't a giant model by 2026 standards. Concrete estimates:

### 8.1 100M_pre at 32 kb context

- **Sequence length post-conv-tower:** 32,768 / 128 = 256 tokens.
- **Attention cost:** O(256² × hidden_dim × num_layers) = O(256² × 768 × 6) ≈ negligible.
- **Conv tower cost:** O(L × hidden_dim) at each level — this dominates. Roughly 5–10× the transformer cost in FLOPs.
- **Memory at batch 8, fp32:** ~6–8 GB activations + 0.4 GB params + optimizer state. Comfortable on a 24 GB GPU.
- **Memory at batch 8, bf16:** ~3–4 GB activations. Comfortable on a 16 GB GPU.
- **Training step time on H100, batch 8, bf16:** roughly 200–400 ms. ~5 step/s. 10k steps in ~30 min. Full epoch over ~10k genomic windows × 50 cells: a few hours.

### 8.2 650M_pre at 32 kb context

- ~5× the param count, ~3× the FLOPs.
- **Memory at batch 4, bf16:** ~12–16 GB activations + 1.4 GB params. Fits on a 24 GB GPU with gradient accumulation.
- **Training step time on H100, batch 4, bf16:** ~700 ms–1 s.
- A full ablation cycle in 1–2 days.

### 8.3 Context length sensitivity

- 32 kb → 131 kb: ~4× memory + 4× time. Still single-GPU feasible at 100M.
- 131 kb → 1 Mb: another 8×. Single-GPU only with gradient checkpointing.

### 8.4 Recommended GPU sizing

- **Smoke test (8M model, 32 kb, batch 32):** any GPU with 16 GB VRAM. Even an RTX 3090 works.
- **Working prototype (100M, 32 kb, batch 16):** 24 GB GPU. A100 40 GB ideal.
- **Scale (650M, 131 kb):** A100 80 GB or H100 80 GB.
- **Full ablation grid:** 4× H100 80 GB is comfortable.

### 8.5 LoRA on top of frozen NTv3

If you eventually fine-tune frozen NTv3 + LoRA (as scooby did on Borzoi):
- LoRA rank 8 on attention QKV + FFN projections: ~0.5–1% of base params trainable.
- **Memory:** essentially the same as frozen forward + small grads. Fits on a 16 GB GPU even for 650M.
- **Speed:** approximately matches frozen forward speed (~1.2× slower with backward).

---

## 9. Existing Work to Crib From

### 9.1 Direct precedents (RNA tracks as input or output, DNA LM backbone)

| Paper / repo | What it does | What you can crib |
|---|---|---|
| **scooby** (Schubach et al. 2024, biorxiv:2024.09.19.613754) — https://github.com/gagneurlab/scooby | Borzoi + LoRA + cell-specific decoder head. Predicts scRNA + scATAC at single-cell resolution from DNA sequence. | The full LoRA recipe for a U-Net-style sequence model. Cell-embedding-based decoder. 8× A40 × 2 days training. |
| **PEFT for Borzoi** (Cofer & Mostafavi 2025, biorxiv:2025.05.26.656171) | Compares LoRA / IA3 / Houlsby on Borzoi for tissue transfer | Concrete PEFT comparison numbers; recommends starting with LoRA rank 8. |
| **Central Dogma Transformer** (arxiv:2601.01089) | Cross-attention from RNA gene representations (Q) to DNA positional representations (K, V) | The cross-attention direction (DNA → RNA, not RNA → DNA). Useful for v1.0 (approach D). |
| **scGPT, Geneformer** | Single-cell foundation models (no DNA) | Cell-embedding patterns; perturbation-prediction objectives. |
| **Enformer / Borzoi** | DNA → multi-track expression prediction | Conv front-end at 128 bp; multi-track loss design. |
| **Personal genome expression** (biorxiv:2025.07.09.664024) | Benchmarks Evo2, NT, Caduceus, Enformer, Borzoi, AlphaGenome on individual-level expression prediction | Sobering negative-results context: pure DNA LMs (Evo2, NT) generally *lose* to Enformer/Borzoi at personal expression. **This is a key motivator for your project.** |

### 9.2 Why scooby is your most direct prior art

scooby:
- Uses a U-Net-style genomic model (Borzoi) — structurally identical to NTv3.
- Uses LoRA — the PEFT method that should transfer cleanly.
- Predicts single-cell coverage at base-pair resolution.
- Provides a cell-embedding-conditioned decoder.

You can essentially "scooby-fy" NTv3:
1. Replace Borzoi → NTv3_100M_pre.
2. Use scooby's cell-embedding generator (Poisson-MultiVI or scPoli) on your Smart-seq data.
3. Apply LoRA to NTv3's transformer interior + add a coverage-prediction head conditioned on cell embedding.

**Caveat re: scooby:** it uses RNA as an *output* target, not an *input* condition. If you want RNA-as-input (which is what the survey doc and your CLAUDE.md describe), scooby gives you the architecture but you'd flip the information flow.

### 9.3 Repos that aren't quite right but worth a look

- `instadeepai/nucleotide-transformer` — the official NTv3 inference notebooks.
- `calico/borzoi` — original Borzoi codebase.
- `gagneurlab/scooby` — *the* most relevant prior art.
- `arc-institute/evo2` — for comparison; Evo2's StripedHyena is a different beast but similar conditioning patterns may apply.

---

## 10. Concrete 1–2 Week Plan

### Day 0 — Environment setup (½ day)

```bash
conda create -n ntv3 python=3.11 -y
conda activate ntv3
pip install "transformers>=4.55.0" "torch>=2.4" accelerate datasets peft scanpy
pip install pyBigWig pysam  # for coverage extraction
huggingface-cli login        # accept InstaDeep license
```

Smoke-test load:

```python
from transformers import AutoTokenizer, AutoModelForMaskedLM
tok = AutoTokenizer.from_pretrained("InstaDeepAI/NTv3_8M_pre", trust_remote_code=True)
m   = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/NTv3_8M_pre", trust_remote_code=True)
out = m(**tok("A"*128, add_special_tokens=False, return_tensors="pt"))
assert out.logits.shape == (1, 128, 11)
```

If this works, the gating / token / dependency story is solved.

### Day 1 — Dataset: BigWig coverage extraction (1 day)

Goal: for a single genotype, produce (per-cell, per-window) coverage at 1 bp and 128 bp resolution.

**Two data sources to consider:**
1. **Jerber et al. (10x, 215 donors)** — *not* full-length but you have it. Use as a placeholder dataset for pipeline development. The pipeline transfers to Smart-seq later.
2. **Cuomo et al. E-MTAB-6945 (Smart-seq2, 125 donors, iPSC → cardiomyocyte)** — the *right* dataset. Per `4_smartseq_datasets_with_genotype_variation.md`. ~36k cells, paired WGS.

For prototype, recommend Cuomo. You need:
- Aligned BAM per cell (or per pool — Smart-seq plate cells are individually barcoded, so re-aligning to per-cell BigWigs is feasible).
- A list of ~1000 32 kb windows tiled across a few chromosomes.
- For each (window, cell): extract bp-resolution coverage, bin to 128 bp, log1p.

Output: `coverage.zarr` or `coverage.npz` of shape `(n_windows, n_cells, 256)` for 32 kb windows binned at 128 bp.

### Day 2 — Implement Approach B (additive input embedding) (1 day)

Code from §4.B above. Add to a script:

```python
# train.py
class NTv3WithRNAInput(nn.Module): ...  # as above

model = NTv3WithRNAInput("InstaDeepAI/NTv3_100M_pre", rna_channels=1)

# Smoke test
ids = torch.randint(5, 9, (2, 32768))           # 32 kb
rna = torch.randn(2, 32768, 1)
out = model(ids, rna)
assert out.logits.shape == (2, 32768, 11)
```

Plus a smoke-test that loss matches vanilla NTv3 when `rna_coverage=0` (zero-init invariant).

### Day 3 — Dataloader + training loop (1 day)

```python
from torch.utils.data import Dataset, DataLoader

class WindowedCoverageDataset(Dataset):
    def __init__(self, windows, ref_genome, coverage_zarr, tokenizer):
        self.windows = windows; self.cov = coverage_zarr; self.tok = tokenizer
        self.ref = ref_genome   # pyfaidx or similar

    def __getitem__(self, i):
        w = self.windows[i]
        seq = self.ref.fetch(w.chrom, w.start, w.end).upper()  # 32k bp string
        ids = self.tok(seq, add_special_tokens=False, return_tensors="pt")["input_ids"][0]
        cells = self.cov[i]  # (n_cells, 256)
        # Pick one cell as the conditioning signal; later we'll batch over cells
        c = cells[np.random.randint(len(cells))]
        # Upsample 256 -> 32k via repeat
        rna_1bp = np.repeat(c, 128)  # (32k,)
        return ids, torch.tensor(rna_1bp[:, None], dtype=torch.float)

# Loss: MLM + (optional) auxiliary coverage prediction
def mlm_loss(model, ids, rna):
    labels, masked_ids = mask_15pct(ids)  # standard 15% mask
    out = model(input_ids=masked_ids, rna_coverage=rna, labels=labels)
    return out.loss
```

Training loop: AdamW, lr=1e-4, weight_decay=0.01, 10k–50k steps. Save checkpoints every 1k steps.

### Day 4 — Smoke runs (1 day)

- Run 1k steps on 8M_pre, 32 kb windows, batch 8. Should fit on any 24 GB GPU.
- Verify: loss decreases. With RNA conditioning zero-init, day-0 loss matches vanilla.
- Verify: `cond.set(zero_tensor)` recovers vanilla loss after training.

### Day 5 — First real experiment (1 day)

Train 100M_pre, 32 kb, batch 16, bf16, on the Cuomo windows × cells dataset. ~24 hours for ~50k steps on a single H100.

Hold out a small set of windows.

### Day 6 — Eval (1 day)

- Compare held-out MLM loss with and without RNA conditioning.
- Compare with shuffled-RNA baseline.
- Compare with right-genotype vs wrong-genotype conditioning. This is the key test.

### Day 7 — Document and decide (½ day)

Write a short results MD. Three outcomes:
1. **Conditioning improves held-out MLM** *and* **right-genotype beats wrong-genotype** → success. Move to approach C and to objective B (genotype-conditioned expression prediction).
2. **Conditioning improves MLM but not genotype-specifically** → information is flowing, but it's a coarse signal. Try approach C for stronger conditioning.
3. **Conditioning doesn't improve** → debug. Likely: RNA is being ignored. Check gradient norms on `rna_proj`. Increase rank of conditioning module.

### Days 8–14 — Approach C + better objective

- Implement FiLM via adaptive layer norm wraps (~3 days).
- Implement Poisson-NLL coverage prediction head (~1 day).
- Train + ablate (~3 days).

End-of-week-2 deliverable: a 100M NTv3 fine-tuned with RNA-conditioned MLM + coverage prediction, evaluated on held-out windows across genotypes, with clear numbers vs the unconditional baseline.

---

## 11. Pitfalls You Will Hit

1. **`inputs_embeds` may not be supported.** Fall back to forward-pre-hooks on the input `Embedding`.
2. **`attention_mask` is implicitly downsampled.** Avoid variable-length batches; pad everything to the fixed window length.
3. **NTv3 is fp16-unsafe in some places** (the conv tower can over/underflow). Train in bf16 if possible; if not, fp32. Check the model card / config for the recommended `dtype`.
4. **Smart-seq coverage is sparse.** Many positions are zero in any one cell. log1p helps but consider adding a "log1p + sign(presence)" two-channel coverage.
5. **Genotype = sequence variation.** If you condition on RNA from individual X, you should also pass X's *sequence* (with its SNPs). Use a personalized reference (e.g., construct via bcftools consensus from VCF). This is the *point* of the project per CLAUDE.md.
6. **Cell axis is variable.** A batch of windows may have different numbers of cells per window. Solve by sampling one cell per window per step (cheapest) or padding (most general).
7. **Catastrophic forgetting.** If you fine-tune NTv3 too aggressively, the MLM ability degrades. Use a small LR (1e-5 for unfrozen, 1e-4 for new params), use LoRA on the existing weights, or simply freeze NTv3 entirely and only train the conditioning module + head.
8. **Mass mismatch.** Don't confuse coverage (read pile-up depth) with TPM (normalized abundance). For Smart-seq full-length, coverage is what you want, normalized only for sequencing depth per cell (CPM-style).
9. **The transformer interior operates at 128 bp resolution.** If you naively inject 1-bp RNA at the transformer layers, the math is wrong. Always match resolutions explicitly.

---

## 12. When To Bail On NTv3 And Pick A Different Backbone

You might reasonably switch to Borzoi or Evo2 if:

- **NTv3's MLM-only pretraining doesn't carry the regulatory grammar.** Borzoi is *supervised* on RNA-seq, so its embeddings have an inductive bias toward regulatory function. The personal-gene-expression benchmark (biorxiv:2025.07.09.664024) found that pure DNA LMs underperform Borzoi on individual-level expression — this is real evidence that supervised models may be stronger for your task.
- **You need 1 Mb context cheap.** Evo2's StripedHyena is faster at long context. NTv3's conv stack helps but isn't free.
- **You want generative capability.** NTv3_generative is already there but the diffusion objective is a specific tool.

For *this prototype*, NTv3 is correct because:
- It has a clean HF interface (Borzoi does not — Borzoi is its own framework).
- The pretrained MLM objective lets you continue training without a bolted-on head.
- The U-Net architecture aligns perfectly with 128 bp coverage binning.

**If after week 2 you have negative results**, my next move would be to *also* run the same recipe on Borzoi (with scooby's LoRA setup) and compare. If Borzoi wins, switch.

---

## 13. Open Questions

1. **What does NTv3's pretrained MLM actually predict well?** The model card touts post-training but is vague on what the pretrained checkpoint alone does. Worth a probing experiment: run NTv3_100M_pre on a held-out human gene body and measure per-position MLM accuracy in exons vs introns vs intergenic. If it's already excellent, your RNA conditioning has less marginal value.
2. **Is the 1 bp embedding stream (16-dim) too narrow a bottleneck for additive conditioning?** It's the *input* embed dim before the stem expands to 768. The conv stem will mix RNA into the model dim, but a 16-dim shift is constrained. May want to inject *after* the stem (at hidden_dim=768) — requires a hook on the stem rather than the embedding.
3. **Can NTv3's species adaptive-LN parameters be "re-trained" for RNA conditioning?** If they're already initialized in the post checkpoint, you could initialize from post and gradient-descend them toward an RNA-driven condition. Cheaper than learning FiLM from scratch.
4. **Does NTv3 do well on the personal-gene-expression benchmark?** The biorxiv 2025.07.09.664024 paper was prior to NTv3's release. Re-running it on NTv3 would give you a baseline to beat.
5. **Should you also evaluate splicing prediction?** Smart-seq covers introns/exons; NTv3 has spliceability per the post checkpoint. A splicing-prediction probe could validate that conditioning carries splicing-relevant signal.

---

## 14. Pointers, Quotes, and Sources

### Model cards (quoted)

- `NTv3_100M_pre`: "U-Net style conv tower → Transformer stack → deconv tower → LM head". "Vocab size 11, Token embedding dim 16, Model (hidden) dim 768, FFN dim 3072, Attention heads 12, Transformer layers 6, Downsample stages 7." "input sequence length need to be a multiple of 128." "Up to 1 Mb of context at nucleotide resolution."
- `NTv3_650M_pre`: "Vocab size 11, Token embedding dim 16, Model (hidden) dim 1536, FFN dim 6144, Attention heads 24, Transformer layers 12, Downsample stages 7."
- `NTv3_8M_pre`: "Vocab size 11, Token embedding dim 16, Model (hidden) dim 256, FFN dim 1024, Attention heads 8, Transformer layers 2, Downsample stages 7."
- `NTv3_100M_post`: "This is a post-trained (supervised) multi-species model that can predict functional genomics tracks (BigWig) and genomic elements (BED) across multiple species."

### Code references

- `nucleotide_transformer_v3/layers.py`: `RotaryEmbedding` class with `UPPER_FREQ = 10000` and Yarn-style rescaling; `MultiHeadAttention` with optional rotary, ESM-style `add_bias_kv`; `SelfAttentionBlock` with pre-LN and GLU FFN; `ConvTowerBlock` (downsample) and `DeconvTowerBlock` (upsample).
- `nucleotide_transformer_v3/adaptive_layers.py`: `AdaptiveLayerNorm` ("Rescale a layer norm with one or several conditions"), `AdaptiveSelfAttentionBlock`, `ConditionedConvTowerBlock`, `ConditionedDeConvTowerBlock`. Zero-init of condition-projection layers.
- `nucleotide_transformer_v3/model.py`: forward returns `{"embedding": ..., "after_transformer_embedding": ..., "logits": ...}`. Inputs: `tokens`, optional `conditions`, optional `species_tokens`.

### Sources

- [InstaDeep NTv3 collection (HuggingFace)](https://huggingface.co/collections/InstaDeepAI/nucleotide-transformer-v3)
- [NTv3 paper (bioRxiv 2025.12.22.695963)](https://www.biorxiv.org/content/10.64898/2025.12.22.695963v1)
- [NTv3 PDF (InstaDeep)](https://instadeep.com/wp-content/uploads/2025/12/NT_v3.pdf)
- [InstaDeep blog: Modelling the Genome with NTv3](https://instadeep.com/2026/02/modelling-the-genome-with-ntv3/)
- [`InstaDeepAI/NTv3_100M_pre` model card](https://huggingface.co/InstaDeepAI/NTv3_100M_pre)
- [`InstaDeepAI/NTv3_650M_pre` model card](https://huggingface.co/InstaDeepAI/NTv3_650M_pre)
- [`InstaDeepAI/NTv3_8M_pre` model card](https://huggingface.co/InstaDeepAI/NTv3_8M_pre)
- [`InstaDeepAI/NTv3_100M_post` model card](https://huggingface.co/InstaDeepAI/NTv3_100M_post)
- [`InstaDeepAI/NTv3_650M_post_131kb` model card](https://huggingface.co/InstaDeepAI/NTv3_650M_post_131kb)
- [InstaDeepAI nucleotide-transformer GitHub](https://github.com/instadeepai/nucleotide-transformer)
- [`nucleotide_transformer_v3.md` docs](https://github.com/instadeepai/nucleotide-transformer/blob/main/docs/nucleotide_transformer_v3.md)
- [DeepWiki NTv3 page](https://deepwiki.com/instadeepai/nucleotide-transformer/3.2-nucleotide-transformer-v3)
- [scooby paper (PMC12615262)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12615262/)
- [scooby code](https://github.com/gagneurlab/scooby)
- [Cofer & Mostafavi 2025: PEFT for Borzoi (bioRxiv 2025.05.26.656171)](https://www.biorxiv.org/content/10.1101/2025.05.26.656171v1)
- [Assessing genomic LMs in personal gene expression (bioRxiv 2025.07.09.664024)](https://www.biorxiv.org/content/10.1101/2025.07.09.664024v1)
- [NVIDIA BioNeMo LoRA blog (Evo2 LoRA target modules)](https://developer.nvidia.com/blog/fine-tuning-biological-foundation-models-with-lora-using-nvidia-bionemo-recipes)
- [Central Dogma Transformer (arXiv 2601.01089)](https://arxiv.org/pdf/2601.01089)

---

## 15. Summary Recommendation

**Start here, in this order:**

1. **Today / day 0:** Accept the InstaDeep license on HF, install transformers >= 4.55, do the smoke load of `NTv3_8M_pre`. Verify `trust_remote_code=True` works in your environment.
2. **Day 1–4:** Implement approach B (additive input embedding with zero-init), build a Cuomo-style coverage dataloader, run MLM-with-RNA on 100M_pre at 32 kb. Validate on held-out windows.
3. **Day 5–7:** Real eval — right-genotype vs wrong-genotype conditioning. Decide whether to scale up.
4. **Day 8–14:** Implement approach C (FiLM via adaptive-LN wraps) for stronger conditioning. Add coverage-prediction head with Poisson NLL.
5. **Week 3+:** Scale to 650M + 131 kb, add cross-attention adapters (approach D), benchmark against scooby and against personal-expression baselines.

**If you only do one thing**: Write the `attach_rna_conditioner` hook from §4.B today. It's <50 lines of code and unblocks every subsequent experiment.

**If you have to make one call against the user's plan:** the MSA-axial-attention idea from `5_architectures_for_rna_conditioned_dna_lms.md` Tier 1 #2 is the *right* architecture for learning gene covariation — but it's a 4–6 week build. The Flamingo-style approach D is the right v1.0. Approach B is the right v0.1. Don't skip B even though it's "just" an input bias.
