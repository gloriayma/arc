# Integrating the Axial-Attention scRNA Architecture with Pretrained NTv3

**Audience:** Gloria Ma (gloria.ma@arcinstitute.org), Arc Institute
**Date:** 2026-06-22
**Scope:** How to slot the axial-attention-over-scRNA-stack architecture (designed in `doing_mds/6.1_variable_depth_strategies.md`) onto a pretrained NTv3 checkpoint. This is a *tactical extension* of `6_quickest_ntv3_prototype_path.md` — that document covers single-cell additive/FiLM conditioning (v0.1/v0.2 in its ladder); this document covers the multi-cell axial-attention upgrade (v0.2/v1.0 in the same ladder).

---

## 1. TL;DR

1. **NTv3 plays the role of the "heavy genomic trunk"** in `6.1`'s recommended architecture. The cell-mixing module (per-cell encoder + row/column axial attention + cell-axis aggregation) sits *before* NTv3's transformer stack and feeds it a single consensus track.
2. **The natural injection point is post-NTv3-conv-tower, at 128 bp resolution.** Bin scRNA coverage to 128 bp from the start. The cell-mixing module operates at 128 bp throughout. Its consensus-track output shares NTv3's transformer-interior resolution 1:1 — no interpolation, no upsampling.
3. **Reuse the FiLM machinery in `adaptive_layers.py`** for the actual injection. NTv3 was already designed to accept FiLM conditioning; the cell-mixing module just provides a richer conditioning signal than a single-cell track.
4. **Build a fresh per-cell encoder** (1-D conv stack with matching 7× downsample). Do not reuse NTv3's conv tower — it expects DNA tokens, not scalar coverage.
5. **Inject at one place initially** (after conv tower, before transformer block 1). Add per-layer FiLM only if signal is too weak. Use zero-init projection so day-0 behavior matches vanilla NTv3.
6. **Pad+mask the cell axis to D=100** per `6.1`. Cell-mixing blocks ignore masked rows in attention softmax; aggregation ignores masked rows.
7. **Stage training:** Phase 1 frozen NTv3 + cell-mixing from scratch → Phase 2 LoRA on NTv3 + cell-mixing full → Phase 3 (only if needed) low-LR full fine-tune.

---

## 2. How This Architecture Maps onto NTv3

### 2.1 The architectural role assignment

From `6.1`, the recommended block-level architecture is:

```
(0) per-cell encoder
(1) cell-mixing blocks × N1=4 (row attention + full column attention)
(2) cell-axis aggregation (attention-pool over D → consensus track)
(3) HEAVY GENOMIC TRUNK ← this is the expensive transformer that operates on a single track
(4) heads
```

**NTv3 _is_ the heavy genomic trunk.** Specifically:
- NTv3 conv tower = pre-trunk downsampling that exists in both architectures.
- NTv3 transformer stack = the trunk in `6.1` step (3).
- NTv3 deconv tower + LM head = the head in `6.1` step (4), already pretrained.

So integration reduces to: (a) building stages (0)–(2) — the cell-mixing preprocessor — and (b) injecting its output into NTv3's transformer stack at the right resolution.

### 2.2 Resolution alignment (the load-bearing constraint)

| Stage | Resolution | Shape |
|---|---|---|
| scRNA coverage input | 1 bp | (B, D, L_bp) |
| After log1p + 128 bp binning | 128 bp | (B, D, L/128, 1) |
| After per-cell encoder | 128 bp | (B, D, L/128, d_cell) |
| After cell-mixing blocks | 128 bp | (B, D, L/128, d_cell) |
| After cell-axis aggregation | 128 bp | (B, L/128, d_cell) |
| After projection d_cell → d_ntv3 | 128 bp | (B, L/128, d_ntv3) |
| NTv3 conv tower output | 128 bp | (B, L/128, d_ntv3) |
| ↑ **THESE TWO MATCH 1:1** ↑ | | |
| NTv3 transformer interior | 128 bp | (B, L/128, d_ntv3) |
| NTv3 deconv output | 1 bp | (B, L, d_ntv3) |
| NTv3 LM head | 1 bp | (B, L, 11) |

The post-cell-mixing consensus track and the post-NTv3-conv-tower hidden state share both resolution (128 bp) *and* position index (token i ↔ genomic bin i). This is the architectural payoff for binning to 128 bp at the start.

For L = 32 kb prototype window: L/128 = 256 tokens — trivial attention cost on both axes (sequence and cell).

### 2.3 What you do NOT do

- ❌ Do not reuse NTv3's conv tower as the per-cell encoder. It expects (V=11, D=16) token inputs. scRNA is a scalar coverage value. Wrong input shape, wrong inductive bias.
- ❌ Do not run NTv3 D=100 times (once per cell). That defeats the entire point of the cell-mixing module and is 100× the compute.
- ❌ Do not feed multi-cell coverage at NTv3's 1 bp input. The conv tower would just average everything together — you'd waste the cell-axis structure. Inject at 128 bp where NTv3's attention can actually use the signal.
- ❌ Do not inject between the conv tower stages. Wait until post-conv-tower (single 128 bp resolution).

---

## 3. Concrete Architecture

### 3.1 Step-by-step walkthrough

The architecture has five logical stages. The first three are new (the cell-mixing preprocessor); the last two are NTv3.

**Stage 0 — Per-cell encoder.** Turn each cell's scalar bp coverage into a feature vector at each 128 bp bin. The same conv weights are applied independently to every cell, with `D` acting as a batch dimension.
- In: `(B, D=100, L_bp, 1)` — one coverage scalar per bp per cell, after `log1p(mean over 128 bp)`.
- Out: `(B, D=100, T, d_cell)` where `T = L_bp/128` and `d_cell ≈ 256`.

**Stage 1 — Cell-mixing blocks (N1 ≈ 4 stacked).** Let cells share information at the same genomic position, and let each cell's track flow information along the genome. Each block applies in order:
- **Row attention along T** (genome axis): reshape `(B, D, T, d) → (B*D, T, d)`, self-attention with RoPE, reshape back. Each cell's track sees itself along the genome.
- **Column attention along D** (cell axis): reshape `(B, D, T, d) → (B*T, D, d)`, self-attention with **no positional encoding** (permutation-invariant), masked padded rows. At each genomic position, all 100 cells attend to each other.
- LayerNorm + MLP per (cell, position) vector.
- In: `(B, D, T, d_cell)` → Out: `(B, D, T, d_cell)` — same shape, enriched.

**Stage 2 — Cell-axis aggregation (the step you asked about).** Collapse the `D=100` axis into a single "consensus" feature vector per genomic position. After this stage the cell axis is gone.
- A small number of learned query vectors (e.g. 1) cross-attend over the `D` cells **independently at each genomic position**:
  - Reshape `(B, D, T, d) → (B*T, D, d)`.
  - For each of the `B*T` (sequence × position) slices, the query attends to the `D` cell vectors at that position, with key-padding-mask from the cell mask.
  - Output one vector per position; reshape back to `(B, T, d_cell)`.
- Mental model: at each position, the model learns "which cells matter most here and how to weight them" — a soft, content-dependent, mask-respecting pool over cells.
- In: `(B, 100, 256, 256)` → Out: `(B, 256, 256)` for the 32 kb prototype.

**Stage 3 — Inject into NTv3 (the "concat" in the abstract `6.1` arch).** Now you have two tensors at the **same** 128 bp resolution and **same** position index:
- `x_dna = (B, T, d_ntv3=768)` — NTv3 hidden state right after the conv tower
- `x_consensus = (B, T, d_cell=256)` — output of Stage 2

These need to be combined and fed into NTv3's transformer stack. `6.1` wrote this loosely as `concat(x_dna, x_consensus)`. In practice, since NTv3 is pretrained at a fixed `d_ntv3 = 768`, **don't actually concat** — that would change input shape and break the pretrained weights. Instead use **additive bias** with a zero-init projection (Option 1 in §4):

```python
x_consensus_projected = Linear(d_cell=256 → d_ntv3=768)(x_consensus)   # zero-init
x_combined = x_dna + x_consensus_projected                              # (B, T, 768)
```

Day-0 model = vanilla NTv3 (because projection is zero). As training proceeds, the projection learns to inject scRNA signal additively. NTv3's transformer stack then processes `x_combined` exactly as it would process `x_dna` alone.

**Stage 4 — NTv3 transformer + deconv + heads.** Run unchanged on `x_combined`. Outputs MLM logits over the input DNA at 1 bp resolution after the deconv tower.

### 3.2 Module diagram

```
DNA: input_ids (B, L_bp)                    scRNA: coverage (B, D=100, L_bp), cell_mask (B, D=100)
    |                                                       |
    v                                            log1p(mean), bin to 128 bp
[NTv3 Embedding (V=11, D=16)]                              |
    |                                                       v
    v                                              (B, D=100, L/128, 1)
[NTv3 Stem conv k=15]                                       |
    |                                                       v
    v                                          [Per-Cell Encoder: 7× downsampling conv stack]
[NTv3 ConvTower: 7 × downsample]                            |          (shared weights across D)
    |                                                       v
    v                                              (B, D=100, L/128, d_cell)
(B, L/128, d_ntv3)                                          |
                                                            v
                                       [Cell-Mixing Block × N1=4]
                                              row attention (L/128 axis, with RoPE)
                                              column attention (D axis, no PE, masked)
                                              LayerNorm, MLP
                                                            |
                                                            v
                                              (B, D=100, L/128, d_cell)
                                                            |
                                                            v
                                          [Cell-Axis Attention Pool over D, masked]
                                                            |
                                                            v
                                              (B, L/128, d_cell)
                                                            |
                                                            v
                                          [Linear d_cell → d_ntv3, zero-init]
                                                            |
                                                            v
                                              consensus_ntv3 (B, L/128, d_ntv3)
    |                                                       |
    |               *** INJECT HERE ***                     |
    +<------------------ FiLM / additive bias ------------- +
    |
    v
[NTv3 Transformer × 6 layers] (conditioned via FiLM or additive bias at each block)
    |
    v
[NTv3 DeconvTower: 7 × upsample, with skip connections from ConvTower]
    |
    v
[NTv3 LM Head] → (B, L, 11) logits for MLM
                + optional auxiliary head → (B, L/128, 1) for coverage prediction
```

### 3.2 Per-cell encoder spec

Goal: scalar coverage at 1 bp → vector at 128 bp, shared weights across cells.

```python
class PerCellEncoder(nn.Module):
    def __init__(self, d_cell=256, num_downsamples=7):
        super().__init__()
        # Match NTv3's downsample count so output resolution matches its transformer interior
        self.stem = nn.Conv1d(1, 32, kernel_size=15, padding=7)
        self.tower = nn.Sequential(*[
            ConvBlock(in_ch, out_ch, downsample=True)
            for in_ch, out_ch in self._channel_schedule(num_downsamples, d_cell)
        ])

    def forward(self, coverage_128bp_upsampled):
        # coverage_128bp_upsampled: (B, D, L_bp, 1)
        # OR pass binned-to-128 directly and skip the tower
        B, D, L, _ = coverage_128bp_upsampled.shape
        x = coverage_128bp_upsampled.view(B * D, L, 1).transpose(1, 2)  # (B*D, 1, L)
        x = self.stem(x)
        x = self.tower(x)  # (B*D, d_cell, L/128)
        return x.transpose(1, 2).view(B, D, L // 128, -1)
```

**Cheaper alternative:** bin coverage to 128 bp first (via `log1p(mean)` over each 128 bp window), then use a much lighter per-cell encoder — just a 2-layer 1-D conv at 128 bp resolution, no downsampling. This is what I'd actually recommend for v0:

```python
class PerCellEncoderLight(nn.Module):
    def __init__(self, d_cell=256, in_channels=1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(in_channels, 128, kernel_size=7, padding=3),
            nn.GELU(),
            nn.Conv1d(128, d_cell, kernel_size=7, padding=3),
        )
    def forward(self, cov_128bp):
        # cov_128bp: (B, D, L/128, C)
        B, D, T, C = cov_128bp.shape
        x = cov_128bp.view(B * D, T, C).transpose(1, 2)
        x = self.net(x)
        return x.transpose(1, 2).view(B, D, T, -1)
```

`d_cell = 256` is a reasonable default — cheaper than `d_ntv3 = 768` and the cell-mixing blocks dominate the param count of this module.

### 3.3 Cell-mixing block spec

```python
class CellMixingBlock(nn.Module):
    def __init__(self, d_cell, n_heads=8):
        super().__init__()
        self.row_attn = MultiHeadSelfAttention(d_cell, n_heads, use_rope=True)
        self.col_attn = MultiHeadSelfAttention(d_cell, n_heads, use_rope=False)  # no PE on cell axis
        self.mlp = nn.Sequential(
            nn.LayerNorm(d_cell), nn.Linear(d_cell, 4 * d_cell),
            nn.GELU(), nn.Linear(4 * d_cell, d_cell),
        )
        self.ln1 = nn.LayerNorm(d_cell); self.ln2 = nn.LayerNorm(d_cell)

    def forward(self, x, cell_mask):
        # x: (B, D, L/128, d_cell); cell_mask: (B, D)  True = real cell, False = pad
        B, D, T, d = x.shape

        # Row attention along L/128, applied per cell (D acts as batch)
        h = x.view(B * D, T, d)
        h = h + self.row_attn(self.ln1(h))                              # (B*D, T, d)
        x = h.view(B, D, T, d)

        # Column attention along D, applied per position (T acts as batch)
        x_col = x.permute(0, 2, 1, 3).contiguous().view(B * T, D, d)    # (B*T, D, d)
        mask_col = cell_mask.unsqueeze(1).expand(B, T, D).reshape(B * T, D)
        x_col = x_col + self.col_attn(self.ln2(x_col), key_padding_mask=~mask_col)
        x = x_col.view(B, T, D, d).permute(0, 2, 1, 3).contiguous()    # (B, D, T, d)

        # MLP
        x = x + self.mlp(x)
        return x
```

Stack `N1 = 4` of these.

### 3.4 Cell-axis aggregation

#### 3.4.1 The problem

After Stage 1, we have a tensor of shape `(B, D=100, T, d_cell)`. We want a single consensus vector per genomic position — shape `(B, T, d_cell)`. So we need to collapse the `D` axis.

How should we summarize 100 cell vectors at a single position into one vector?

| Option | Pros | Cons |
|---|---|---|
| **Mean pool** `out = mean_i(v_i)` | Trivial, permutation-invariant, depth-agnostic | Every cell weighted equally regardless of content. Noisy cells contribute as much as informative ones. Throws away covariance. |
| **Max pool** `out_j = max_i(v_{ij})` | Picks the most extreme value per dim | Loses everything but the extreme. Outlier-sensitive. |
| **Sum + MLP (DeepSets)** | Universal approximator for permutation-invariant functions | In practice less sample-efficient than attention pool. Good side-channel only. |
| **Attention pool with learned query** | Model **learns which cells matter** at each position. Permutation-invariant. Same weights at any `D`. | Slightly more parameters. |
| **Concat all D vectors → MLP** | Maximally expressive | **Breaks permutation invariance.** Forces fixed `D`. Hard rejection — violates load-bearing property from `6.1`. |

We use **attention pool with learned queries**. The next subsection explains exactly what that means.

#### 3.4.2 What "learned query vectors" actually means

Standard cross-attention has three inputs: queries `Q`, keys `K`, values `V`. The operation is:
```
attention(Q, K, V) = softmax(Q · Kᵀ / √d) · V
```
"Queries" are what we're asking; "keys" determine relevance; "values" are what we extract.

In **normal** cross-attention (e.g., decoder attending to encoder), `Q` comes from one sequence and `K, V` from another — both are derived from data.

In **learned-query** cross-attention, `Q` is a **parameter of the model** — a `(k, d)` tensor of learnable vectors. There are `k` of them (often `k=1`). The same `Q` is used for every batch element. The model learns what to "look for" in the input set.

Concretely for our case at position `t`:

```
Inputs at this position:
  cell vectors  v_1, v_2, ..., v_100 ∈ R^{d_cell}    (from stage 1, padded rows masked)

Learned parameters of the CellAxisPool module:
  q ∈ R^{d_cell}                                       (k=1 query)

Operation:
  scores_i = (q · v_i) / sqrt(d_cell)                  for i = 1..100
  weights  = softmax(scores)   (over the 100 cells)    ← masked rows get -inf score
  output   = sum_i weights_i · v_i                     ← (d_cell,) summary

  Mental model: q is a "question" the model has learned to ask
  about cells. Cells whose feature vectors resemble q (high
  dot product) get more weight. The summary is a weighted
  combination of value vectors where the weights are content-dependent.
```

This is applied **independently at every genomic position `t`** — the same query parameter `q`, fresh attention weights for the cells at that position.

**Why this preserves permutation invariance on `D`:** the softmax is computed over the `D` axis. Permuting the cell vectors permutes the `scores_i` correspondingly; the `softmax` weight assigned to each cell `i` depends only on `v_i` (not its position in the list); the final sum is unchanged by permutation. This is the same property that makes MSA Transformer's column attention permutation-equivariant — it's a sum-over-keys operation.

**Why one query is often enough:** for a single consensus track, you want one summary vector per position. Multiple queries (`k > 1`) give you `k` summaries; you can average them, concat them, or treat them as parallel "expert" extractors. Multiple queries can capture different aspects of the cell set (e.g., one query learns to extract "mean-like" signal, another extracts "high-variance positions") but in practice the model has to discover this and it doesn't always do so cleanly. Start with `k=1`; ablate up.

#### 3.4.3 Previous art

The "learned query cross-attention" pattern has many independent inventions and many names. The most directly relevant references:

| Paper | What they called it | Connection |
|---|---|---|
| Devlin et al. 2019 (**BERT**) | `[CLS]` token | The original pattern: prepend a learned token to the sequence; its final-layer representation is the sequence summary. Single learned query, applied via self-attention. Conceptually identical to single-query attention pool. |
| Dosovitskiy et al. 2021 (**ViT**) | class token | Same pattern, image patches instead of word tokens. |
| **Lee et al. 2019 (Set Transformer)** | **PMA — Pooling by Multihead Attention** | **The most direct precedent — literally what we use.** Define `k` learnable "seed vectors" `S`; cross-attend `S` as queries to the input set. `PMA_k(Z) = MAB(S, Z)`. Proven permutation-invariant on the input set. Theoretical universal-approximation results for permutation-invariant functions. |
| Carion et al. 2020 (**DETR**) | object queries | 100 learnable queries cross-attend to image features; each query learns to focus on a different object. Output is detection-by-query. |
| Jaegle et al. 2021/2022 (**Perceiver / Perceiver IO**) | latent array | A large learnable latent (e.g., 256–1024 vectors) cross-attends to giant variable-size inputs; decouples downstream compute from input size. |
| Alayrac et al. 2022 (**Flamingo**) | Perceiver Resampler | 64 learned latents bridge variable visual features into a frozen language model via cross-attention. |
| Vaswani et al. 2017 (**original Transformer**) | decoder cross-attention | The encoder-decoder cross-attention is the parent operation; what we use specializes the query side to be learned parameters rather than data-derived. |

**The right citation for us is Set Transformer's PMA.** It is the cleanest framing of "fixed-size summary of a permutation-invariant set via cross-attention with learned seeds," and the permutation-equivariance proof in §3 of the paper is exactly what we need. The Set Transformer paper is also where ISAB comes from (which we considered and rejected in `6.1` because `D ≤ 100 < M=128` — the bottleneck-vs-full-attention crossover doesn't favor ISAB at our scale, but PMA still does the right thing for aggregation).

#### 3.4.4 Implementation

```python
class CellAxisPool(nn.Module):
    """Set Transformer-style PMA over the cell axis at each genomic position."""
    def __init__(self, d_cell, n_queries=1, n_heads=8):
        super().__init__()
        # Learned query vectors — the "seeds" in Set Transformer terminology
        self.queries = nn.Parameter(torch.randn(n_queries, d_cell) * 0.02)
        self.attn = nn.MultiheadAttention(d_cell, n_heads, batch_first=True)

    def forward(self, x, cell_mask):
        # x: (B, D, T, d_cell); cell_mask: (B, D), True = real cell
        B, D, T, d = x.shape

        # Move T into the batch dim so the pool runs independently per position
        x_in = x.permute(0, 2, 1, 3).contiguous().view(B * T, D, d)

        # Broadcast the same learned queries across all (B, T) positions
        q = self.queries.unsqueeze(0).expand(B * T, -1, -1)               # (B*T, n_queries, d_cell)

        # Mask out padded cell rows
        mask = ~cell_mask.unsqueeze(1).expand(B, T, D).reshape(B * T, D)  # True = ignore

        # Cross-attention: queries attend over the D cells at each position
        out, _ = self.attn(q, x_in, x_in, key_padding_mask=mask)          # (B*T, n_queries, d_cell)

        # If n_queries > 1, average them into a single consensus vector per position
        out = out.mean(dim=1).view(B, T, d)                                # (B, T, d_cell)
        return out
```

Output: `(B, L/128, d_cell)` consensus track. This feeds Stage 3 (NTv3 injection).

### 3.5 Projection + injection

```python
self.consensus_proj = nn.Linear(d_cell, d_ntv3)
nn.init.zeros_(self.consensus_proj.weight)  # zero-init: day-0 model = vanilla NTv3
nn.init.zeros_(self.consensus_proj.bias)
```

Inject via one of three mechanisms (see §4).

---

## 4. Injection Mechanisms — Three Options

In order of complexity. All preserve day-0-equals-vanilla-NTv3 with zero-init.

### Option 1 — Additive bias on transformer input (simplest, ~30 LOC)

Add `consensus_ntv3` once to the hidden state at the input of NTv3's transformer stack:

```python
def forward(self, dna_ids, rna_cov, cell_mask, labels=None):
    consensus = self.cell_mixing(rna_cov, cell_mask)         # (B, L/128, d_cell)
    consensus = self.consensus_proj(consensus)               # (B, L/128, d_ntv3), zero-init

    # Run NTv3 conv tower
    h = self.ntv3.embedding(dna_ids)
    h = self.ntv3.stem(h)
    h = self.ntv3.conv_tower(h)                              # (B, L/128, d_ntv3)

    # *** INJECT ***
    h = h + consensus

    # Continue NTv3
    h = self.ntv3.transformer(h)
    h = self.ntv3.deconv_tower(h, conv_skips)
    return self.ntv3.lm_head(h)
```

**Pros:** dead-simple, debuggable, zero-init invariant trivial.
**Cons:** signal injected once; relies on NTv3 transformer to propagate it across all layers.
**Use first.** If MLM improvement is decisive, may not need anything fancier.

### Option 2 — FiLM via NTv3's `AdaptiveLayerNorm` (recommended v0.2, ~100 LOC)

NTv3 ships `adaptive_layers.py` with `AdaptiveLayerNorm`, `AdaptiveSelfAttentionBlock`, `ConditionedConvTowerBlock`, and `ConditionedDeConvTowerBlock`. The posttrained checkpoint uses these for species conditioning. Repurpose for RNA conditioning.

```python
class FiLMFromConsensus(nn.Module):
    """Per-token FiLM: scale and shift driven by consensus track."""
    def __init__(self, d_ntv3):
        super().__init__()
        self.proj = nn.Linear(d_ntv3, 2 * d_ntv3)
        nn.init.zeros_(self.proj.weight); nn.init.zeros_(self.proj.bias)

    def forward(self, h, consensus):
        # h, consensus: (B, L/128, d_ntv3)
        gamma_beta = self.proj(consensus)
        gamma, beta = gamma_beta.chunk(2, dim=-1)
        return h * (1 + gamma) + beta
```

Wrap each transformer-block LN with a FiLM call. Same `consensus` tensor reused at every block.

**Pros:** Conditioning at every block; uses NTv3's existing inductive bias.
**Cons:** More code; need to find the right LN instances to wrap.
**Use:** if Option 1 underfits.

### Option 3 — Gated cross-attention (Flamingo-style, ~200 LOC)

Insert new cross-attention modules between NTv3 transformer blocks:

```python
class GatedCrossAttnAdapter(nn.Module):
    def __init__(self, d_ntv3, n_heads=8):
        super().__init__()
        self.attn = nn.MultiheadAttention(d_ntv3, n_heads, batch_first=True)
        self.gate = nn.Parameter(torch.zeros(1))  # zero-init gate

    def forward(self, h, consensus):
        # h, consensus: (B, L/128, d_ntv3)
        attn_out, _ = self.attn(h, consensus, consensus)
        return h + torch.tanh(self.gate) * attn_out  # tanh(0) = 0
```

**Pros:** Most expressive; each NTv3 transformer layer can independently query the consensus track.
**Cons:** Most code; biggest deviation from pretrained NTv3 forward path.
**Use:** v1.0 after Options 1/2 have been compared.

---

## 5. Where Specifically in NTv3's Code to Hook

Per the existing NTv3 prototype doc, NTv3 exposes `model.core` and supports config-driven hooks. Concrete options:

**Forward pre-hook on the transformer stack** (Option 1, simplest):
```python
def transformer_input_hook(module, args, kwargs):
    h = args[0]
    h = h + self._cached_consensus
    return (h,) + args[1:], kwargs

ntv3.core.transformer.register_forward_pre_hook(
    transformer_input_hook, with_kwargs=True
)
```

**Wrap individual LayerNorm modules** (Option 2):
```python
for name, module in ntv3.core.transformer.named_modules():
    if isinstance(module, nn.LayerNorm) and "norm" in name:
        wrapper = FiLMWrap(module, d_ntv3, ...)
        parent_name, attr = name.rsplit(".", 1)
        parent = ntv3.core.transformer.get_submodule(parent_name)
        setattr(parent, attr, wrapper)
```

**Insert cross-attn between blocks** (Option 3): requires either monkey-patching `transformer.forward` or forking the file. Hooks aren't sufficient here because cross-attn needs to run *between* blocks, not before/after them.

**Verification step before building any of this:** read `ntv3_huggingface_new.py` (loaded via `trust_remote_code=True`), find the transformer stack's `forward()`, confirm `inputs_embeds` keyword and hook patterns work as expected. The existing prototype doc §2.5 flagged some ambiguity here — check before you commit code.

---

## 6. Training Schedule

Modified from the existing prototype plan §10:

### Phase 0 — Smoke test (½ day)
- Build the full module assembled with zero-init projection.
- Forward pass with `rna_coverage` = zeros → output should match vanilla NTv3 *exactly* (bit-for-bit if zero-init is correct).
- This is the critical invariant. If it fails, debug before anything else.

### Phase 1 — Cell-mixing module training with frozen NTv3 (~3–5 days)
- Freeze all NTv3 parameters. Only train: per-cell encoder, cell-mixing blocks, cell-axis pool, consensus projection.
- Loss: MLM (regularizer) + auxiliary masked-coverage prediction on held-out cells.
- Trainable params: ~5–20M (much smaller than NTv3 100M).
- Compute: single H100, batch 16, 32 kb window, bf16. ~50k steps in ~24 hours.
- **Validates:** the cell-mixing module is capable of producing useful conditioning signal, given NTv3's pretrained sequence understanding.

### Phase 2 — LoRA on NTv3 + cell-mixing full (~3–5 days)
- Add LoRA (rank 8) on NTv3's transformer Q/K/V/FFN projections.
- Continue training cell-mixing module + LoRA.
- Allows NTv3 to adapt slightly without catastrophic forgetting.
- Compute: same as Phase 1, slightly slower (~20%).
- **Validates:** that letting NTv3 adapt yields further gain.

### Phase 3 — Optional low-LR full fine-tune (~1 week)
- Unfreeze all of NTv3 at lr=1e-5; keep cell-mixing at lr=1e-4.
- Only do this if Phase 2 plateaus and you have a clear hypothesis about why NTv3 needs to adapt further.
- High risk of catastrophic forgetting — monitor MLM loss carefully.

### Evaluation at each phase
1. **Held-out MLM loss** with vs. without RNA conditioning (zero out `rna_cov` at eval).
2. **Right-genotype vs wrong-genotype conditioning:** for held-out windows, condition on (a) the correct individual's scRNA, (b) a randomly chosen other individual's scRNA. Right-genotype should outperform.
3. **Coverage-prediction held-out cells:** mask one cell at training time, predict its coverage from the others. R² or Pearson on log1p coverage.
4. **Promoter-enhancer pairing benchmarks** (e.g., ENCODE-rE2G, CRISPRi-FlowFISH ground truth): measure attention/attribution between annotated enhancer–TSS pairs vs. random pairs.

---

## 7. Hyperparameters

Starting point — adjust based on smoke-test loss curves.

| Parameter | Value | Reasoning |
|---|---|---|
| NTv3 checkpoint | `InstaDeepAI/NTv3_100M_pre` | Per existing prototype doc; 100M = best size/iteration trade-off |
| `d_ntv3` | 768 | NTv3 100M hidden dim |
| `d_cell` | 256 | Cheap; cell-mixing blocks dominate cell-side params |
| `N1` (cell-mixing blocks) | 4 | Two rounds of cross-cell info flow |
| `D_max` | 100 | Data ceiling per `6.1` |
| `D_actual` distribution | `Uniform[1, 100]` per step | Depth-robust weights |
| `L_bp` (window) | 32,768 (32 kb) | NTv3 native benchmark length |
| `L/128` (tokens) | 256 | Trivial attention on both axes |
| Coverage binning | log1p(mean over 128 bp) per cell | Standard scRNA / Enformer / Borzoi |
| Row PE in cell-mixing | RoPE | Match NTv3 |
| Column PE in cell-mixing | none | Permutation invariance on cell axis |
| Injection mechanism | additive bias (Option 1) | Start simple |
| Consensus projection init | zero | Day-0 = vanilla NTv3 invariant |
| Cell-mixing learning rate | 1e-4 | New params |
| NTv3 learning rate, Phase 1 | 0 (frozen) | |
| NTv3 LoRA learning rate, Phase 2 | 1e-4 | |
| NTv3 full LR, Phase 3 | 1e-5 | Avoid catastrophic forgetting |
| Batch size | 16 (bf16) | H100 80 GB |
| MLM loss weight | 1.0 | NTv3's pretraining loss |
| Coverage-prediction loss weight | 1.0 | Auxiliary objective |
| Optimizer | AdamW, wd=0.01 | NTv3 default |

---

## 8. Decisions Deferred / Open Questions

1. **Inject at one location or every transformer block?** Start with one (after conv tower); upgrade to per-block FiLM only if signal is too weak. Empirical question.
2. **Per-cell metadata conditioning.** Cell type, donor ID — add as additive vectors to the per-cell encoder output, or omit. I lean: add, since each cell does carry meaningful identity, but keep the cell-axis itself permutation-invariant. Worth ablating.
3. **Two-channel coverage** (mean + max, or coverage + junction reads). Smart-seq junction reads carry splicing signal that mean-coverage hides. Borzoi uses junction tracks; consider adding as a second input channel.
4. **Auxiliary head architecture.** Linear from final NTv3 hidden state to 128 bp coverage prediction. Could also be Borzoi-style 2-layer head. Pick simple first.
5. **Should the consensus track also enter NTv3 *before* the conv tower (1 bp resolution)?** §6.B / §6.C in the existing doc suggest this for single-cell. With a multi-cell consensus you could upsample via `repeat_interleave(128)`. Skip for v0; revisit if interior-only injection isn't enough.
6. **Resolution of the 5-downsample variant** (32 bp interior, `NTv3_5downsample_pre`). Higher resolution may help splicing-related signal. Try after the 7-downsample baseline works.
7. **D = 1 fallback behavior.** With `Uniform[1, 100]` training, model should degrade gracefully to DNA-only. Verify with explicit eval.

---

## 9. What This Does NOT Cover

- **Data extraction pipeline** for per-cell coverage from Smart-seq BAM files → `coverage.zarr`. Covered in existing prototype doc §10 Day 1 and in `doing_mds/3.x` alignment pipeline notes.
- **Choice of dataset** (Jerber vs Cuomo). Covered in `doing_mds/1_jerber_data.md` and `doing_mds/2_cuomo_data.md`.
- **Personal genome construction** (VCF → individual reference). Covered in `doing_mds/2.3_diploid_vs_haploid.md` and `doing_mds/2.4_diploid_alignment.md`.
- **Promoter-enhancer eval datasets and benchmarks.** Separate research thread; see existing CURRENT_LIMS-style docs.
- **Alternative backbones** (Borzoi, Evo2). Existing prototype doc §12 covers when to switch.

---

## 10. Relationship to Other Planning Docs

| Doc | Role | Relationship |
|---|---|---|
| `planning_mds/5_architectures_for_rna_conditioned_dna_lms.md` | Survey of architectural space | This doc instantiates the MSA / axial-attention path from Tier 1 #2 of that survey |
| `planning_mds/6_quickest_ntv3_prototype_path.md` | NTv3 fine-tuning recipe (single-cell conditioning) | This doc extends approaches B/C/D from §4 of that doc to multi-cell axial attention |
| `doing_mds/6_axial_genomics_precedents.md` | Literature search on axial attention in genomics | Gap analysis — establishes this is white space |
| `doing_mds/6.1_variable_depth_strategies.md` | Variable-depth strategy comparison | Provides the cell-axis design (pad+mask, full column attention, no PE) at D ≤ 100 |
| `doing_mds/5.1_ntv3_arch.html` | NTv3 architecture reference | Source-of-truth for layer dims, conv-tower depth, FiLM machinery |

---

## 11. End-to-End Code Skeleton

```python
import torch
import torch.nn as nn
from transformers import AutoModelForMaskedLM


class NTv3WithAxialAttention(nn.Module):
    def __init__(self, repo="InstaDeepAI/NTv3_100M_pre", d_cell=256, N1=4):
        super().__init__()
        self.ntv3 = AutoModelForMaskedLM.from_pretrained(repo, trust_remote_code=True)
        d_ntv3 = self.ntv3.config.hidden_size  # 768 for 100M

        self.cell_encoder = PerCellEncoderLight(d_cell=d_cell)
        self.cell_mixing = nn.ModuleList([
            CellMixingBlock(d_cell=d_cell) for _ in range(N1)
        ])
        self.cell_pool = CellAxisPool(d_cell=d_cell)
        self.consensus_proj = nn.Linear(d_cell, d_ntv3)
        nn.init.zeros_(self.consensus_proj.weight)
        nn.init.zeros_(self.consensus_proj.bias)

        self._cached_consensus = None
        self._register_injection_hook()

    def _register_injection_hook(self):
        def pre_hook(module, args, kwargs):
            if self._cached_consensus is not None:
                args = (args[0] + self._cached_consensus,) + args[1:]
            return args, kwargs
        # Replace with actual transformer-stack attribute path after inspecting NTv3:
        self.ntv3.core.transformer.register_forward_pre_hook(pre_hook, with_kwargs=True)

    def encode_cells(self, rna_cov_128, cell_mask):
        # rna_cov_128: (B, D, L/128, 1)  ;  cell_mask: (B, D)
        x = self.cell_encoder(rna_cov_128)             # (B, D, T, d_cell)
        for block in self.cell_mixing:
            x = block(x, cell_mask)
        consensus = self.cell_pool(x, cell_mask)       # (B, T, d_cell)
        return self.consensus_proj(consensus)          # (B, T, d_ntv3)

    def forward(self, input_ids, rna_cov_128, cell_mask, labels=None, **kw):
        self._cached_consensus = self.encode_cells(rna_cov_128, cell_mask)
        try:
            out = self.ntv3(input_ids=input_ids, labels=labels, **kw)
        finally:
            self._cached_consensus = None
        return out
```

This is a single-file skeleton, ~150 LOC excluding the encoder/mixing/pool module definitions. The full set (including those modules) is ~400 LOC.

---

## 12. Critical Path to a First Result

1. **Day 0:** Smoke-load `NTv3_100M_pre`. Inspect `ntv3_huggingface_new.py` to locate the transformer-stack attribute and confirm a forward pre-hook works. (~½ day)
2. **Day 1–2:** Implement the per-cell encoder + cell-mixing blocks + pool + projection. Unit test with random input. Verify zero-init invariant. (~2 days)
3. **Day 3:** Wire up the dataloader (extend the existing prototype doc's `WindowedCoverageDataset` to return D cells per window). (~1 day)
4. **Day 4–6:** Phase 1 training run (frozen NTv3 + cell-mixing) on Cuomo Smart-seq, 32 kb windows. (~3 days including monitoring)
5. **Day 7:** Eval: held-out MLM with/without conditioning + right-genotype vs wrong-genotype. (~1 day)
6. **Day 8+:** Decide based on Phase 1 results whether to (a) move to Phase 2 LoRA, (b) upgrade to FiLM injection, or (c) debug if no signal.

End-of-week-2 deliverable: a working axial-attention NTv3 with first-pass numbers vs. the single-cell conditioning baseline from the existing prototype doc.
