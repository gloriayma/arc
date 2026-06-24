# Incorporating RNA Expression Signals into DNA Language Models

**Context:** The goal is to inject scRNA-seq coverage tracks (BigWig-style, full base-pair resolution, Smart-seq) into DNA LMs like NTv3 and Evo2, with a primary objective of learning **covariation between genes** across genotypes/conditions. This document surveys architectural options.

---

## 1. MSA-Inspired Approaches — The Researcher's Framing

### The analogy

In protein biology, a Multiple Sequence Alignment presents homologous sequences as a matrix of shape `(num_sequences, sequence_length)`. The key insight is that **columns co-vary**: positions that interact physically tend to mutate together across the evolutionary record. ESM-MSA-1b (Rao et al. 2021, "MSA Transformer") exploits this with **axial attention** — alternating attention along rows (sequence positions within one sequence) and along columns (the same position across all sequences in the alignment). This is cheaper than full 2D attention over the flattened matrix (O(M·L) per head instead of O((M·L)²)).

Your setup maps onto this as:
- **Rows** = cells or samples (genotypes, conditions, cell types)
- **Columns** = genomic positions
- **Values** = RNA coverage at that position

This is structurally identical to an MSA, and the co-evolutionary signal you want to learn (covariation of gene expression across genotypes) is exactly the "columnar covariation" the MSA Transformer was designed for.

### Precedents and related work

**ESM-MSA-1b / MSA Transformer (Rao et al. 2021, NeurIPS)**
- Axial attention: row-wise (within a sequence) interleaved with column-wise (across sequences at the same position).
- Column attention is cheap because the sequences are short (proteins, ~few hundred residues). Genomic sequences are 3–6 orders of magnitude longer, which is the core scaling challenge.
- Key result: column attention is what captures coevolutionary constraints. Without it, performance degrades to single-sequence quality.

**AlphaFold2 MSA Stack (Jumper et al. 2021, Nature)**
- More elaborate: MSA row attention with column attention, plus triangular multiplicative updates for pair representation.
- The pair representation is a separate `(L, L)` matrix that gets updated by the MSA. For DNA this would be a `(L, L)` genomic contact map, which maps naturally onto Hi-C data — a potential extension.
- The triangular updates ("triangle self-attention", "triangle multiplicative update") enforce transitivity: if position A contacts B, and B contacts C, that influences A-C. This is directly relevant for enhancer-promoter chains.

**Perceiver / Perceiver IO (Jaegle et al. 2021)**
- A cross-attention mechanism that maps an arbitrarily large input array down to a fixed-size latent array, then processes in latent space.
- Could be used to compress M × L RNA tracks into a fixed latent, then cross-attend from the DNA sequence to the latent.

### Tradeoffs of axial attention for position × sample

| Consideration | Notes |
|---|---|
| Complexity | Row attention: O(M · L²). Column attention: O(L · M²). If L = 500k and M = 100 cells, column attention is negligible (100² = 10k) but row attention is massive. |
| Long-context wall | At full base-pair resolution over a gene body (say L = 200k), standard attention over L is infeasible. Must use sparse attention, sliding windows, or hierarchical pooling. |
| Position alignment | Requires all tracks mapped to the same reference coordinates — standard for BigWig against a reference genome, but care needed for multi-genotype data with structural variants. |
| Positional encoding | Genomic distance encodes regulatory grammar. Relative positional encodings (ALiBi, RoPE, Enformer's convolution front-end) matter more than absolute. |
| Sparse signal | Single-cell coverage is sparse. Most positions will be zero for most cells. Need explicit handling (masking, special zero-token). |

### Practical recommendation on the MSA framing

The MSA framing is sound and has strong precedent. The main adaptation needed is **long-context handling**. A realistic approach:
1. Use a **convolutional front-end** (like Enformer's) to pool coverage into tokens at ~128 bp resolution, reducing L by 128×.
2. Apply axial attention on the resulting `(M, L/128)` matrix.
3. Optionally reinject fine-grained sequence features from the DNA LM at the same resolution.

---

## 2. Cross-Modal Attention — Conditioning a DNA LM on RNA Tracks

### The framing

Instead of treating RNA tracks as extra "rows" of a joint matrix, treat them as a **separate modality** that conditions the DNA LM via cross-attention, analogous to how:
- Text-to-image models (Stable Diffusion, DALL-E 2) condition image generation on text embeddings via cross-attention injected into the U-Net.
- Flamingo (Alayrac et al. 2022) conditions a frozen language LM on visual tokens via gated cross-attention layers.

The RNA tracks become a "context" that the DNA model queries into. The DNA model processes sequence tokens; at each layer (or every K layers), the DNA tokens cross-attend to a compressed RNA representation.

### Architecture sketch

```
RNA coverage tracks (M cells × L positions)
          |
   [RNA Encoder: CNN or Transformer or Hyena/Mamba]
          |
   RNA latent tokens: shape (N_rna_tokens, D)   ← compressed
          |
          | cross-attention
          v
DNA LM hidden states (L_dna tokens, D)
   — standard DNA LM self-attention layers
   — interleaved with cross-attention to RNA latent
```

### Relevant precedents

**Flamingo (Alayrac et al. 2022, DeepMind / NeurIPS)**
- Frozen vision encoder + frozen LLM + learned cross-attention layers inserted between frozen LM layers.
- Demonstrates that a frozen pretrained LM can be conditioned on a new modality with minimal new parameters.
- Directly applicable: freeze Evo2 or NTv3, insert cross-attention layers that attend to RNA embeddings.

**Perceiver IO (Jaegle et al. 2021)**
- Cross-attention from a small learned query array to a large input (RNA tracks), producing fixed-size output regardless of input size. Solves the scalability problem: the RNA encoder compresses M × L into N fixed tokens.

**Hyena / Mamba for the RNA encoder**
- If RNA tracks are very long (L = 500k positions), a transformer encoder is expensive. State space models (Mamba, Hyena) scale linearly in L and could serve as the RNA encoder.
- Evo2 itself uses Hyena convolutions. The same mechanism could be used to encode the RNA side, creating architectural consistency.

### Tradeoffs vs. MSA framing

| | MSA / Axial Attention | Cross-Modal Attention |
|---|---|---|
| Coupling between modalities | Tight: DNA and RNA tokens attend to each other at every layer | Looser: RNA compressed to latent first, then injected |
| Parameter cost | High (full joint model) | Lower if RNA encoder is separate |
| Preserves pretrained DNA LM | Requires significant retraining | Flamingo-style: can freeze DNA LM |
| Fine-grained position-level RNA features | Natural | Requires careful positional encoding in RNA encoder |
| Learning covariation across cells | Yes (column attention) | Depends on whether RNA encoder captures cross-cell variation |

Cross-modal attention is more **modular** and better for leveraging a pretrained DNA LM. The MSA approach is more **principled** for learning positional covariation.

---

## 3. Multi-Task / Expression-Supervised DNA Models

### The paradigm

These models treat expression prediction as a **supervised training signal**: given a DNA sequence, predict the RNA output. They do not take RNA as input — RNA shapes the loss, not the input — but they reveal critical architectural lessons.

### Key models

**Enformer (Avsec et al. 2021, Nature Methods)**
- Input: 196,608 bp of DNA sequence.
- Architecture: convolutional front-end (reduces L from 196k to 1536 tokens at 128 bp resolution), followed by Transformer encoder, followed by pointwise prediction head per track.
- Predicts ~5000 genomic tracks (CAGE, DNase, ChIP-seq, RNA-seq) from 4,000+ cell types/tissues.
- Key insight: the Transformer's self-attention learns to associate distal regulatory elements with their target genes — it implicitly learns enhancer-promoter grammar from the multi-track supervised signal.
- Limitation: predicts population-average expression from the reference genome. Does not model genotype variation.

**Borzoi (Linder et al. 2023, Nature Genetics)**
- Extends Enformer to RNA-seq at single-nucleotide resolution (not just CAGE peaks). Uses a U-Net-style decoder to upsample predictions back to base-pair resolution.
- More directly analogous to your BigWig coverage tracks.
- The U-Net skip connections allow fine-grained position predictions without losing long-range context.
- Scales to 524,288 bp input.

**Sei (Chen et al. 2022, Nature Genetics)**
- 21,907 chromatin profiles as prediction targets. Uses a "sequence class" system to generalize across cell types.

**Orca (Zhou 2022, Nature Methods)**
- Predicts Hi-C chromatin contact maps from sequence. Relevant if incorporating 3D genome structure (chromatin looping is key to enhancer-promoter communication).

### What this paradigm teaches

1. **The supervised signal from multi-track RNA/chromatin data is enormously powerful for learning regulatory grammar.** Enformer's emergent ability to localize enhancers comes from this supervision, not from regulatory annotations.
2. **128 bp tokenization** is the right resolution for most regulatory signals; base-pair resolution is needed only for splice sites and specific TF binding motifs.
3. **Multi-task training over many cell types** forces learning of cell-type-specific regulatory programs — directly relevant to your goal.
4. **Predicting tracks is not the same as conditioning on them.** Reversing the information flow (RNA → sequence understanding, not sequence → RNA prediction) is a different and harder problem.

### The "inverse Enformer" idea

An underexplored direction: train a model to **infer regulatory state from RNA tracks** (the inverse of Enformer) and use those inferred states as conditioning signals for the DNA LM. This is learning a "cell state embedding" from RNA coverage, then using it as context.

---

## 4. Joint Pretraining — DNA + RNA from Scratch vs. Finetuning

### Joint pretraining from scratch

**Precedents in biology:**
- **Geneformer (Theodoris et al. 2023, Nature):** Pretrained on 30M single-cell transcriptomes. Uses gene rank order (by expression) as input tokens. Not DNA-sequence-aware, but demonstrates large-scale scRNA pretraining produces powerful cell-state representations.
- **scGPT (Wang et al. 2023, Nature Methods):** GPT-style pretraining on scRNA-seq with gene tokens and expression value tokens. Enables transfer across cell types and perturbation prediction.
- **UCE (Rosen et al. 2023):** Universal Cell Embedding using ESM2 protein embeddings as proxies for gene identity — an early attempt to bridge sequence and expression.
- **Evo2 (Arc Institute / Nguyen et al. 2024):** 7B parameter model trained on 9.3T nucleotides. Uses StripedHyena (Hyena + attention) at 1M context. Does not incorporate RNA data but is a key base model for finetuning.

**Tradeoffs:**
- **Pro:** The model can learn joint representations from the ground up. No modality mismatch from misaligned pretraining objectives.
- **Con:** Enormously expensive. Requires aligning DNA sequence tokens with RNA coverage tokens at the same genomic coordinates. Data curation is a major bottleneck.
- **Con:** Current DNA LMs (Evo2, NTv3) represent years of compute investment. Discarding this for joint pretraining from scratch is only worthwhile with comparable compute.

### Finetuning a pretrained DNA LM

More practical. Two main variants:

**a) Full finetuning with RNA as auxiliary input**
- Add RNA track tokens to the input (embed coverage values as continuous tokens, concatenate with DNA tokens).
- Continue pretraining or finetune on a task-specific loss.
- Risk: catastrophic forgetting of DNA pretraining knowledge.
- Mitigation: low learning rate, LoRA/adapter layers, replay of DNA pretraining data.

**b) Frozen DNA LM + learned RNA encoder + task head (Flamingo paradigm)**
- Freeze Evo2 / NTv3 weights. Train only: RNA encoder, cross-attention adapter layers, task prediction head.
- Much cheaper. Preserves DNA sequence knowledge. Limits how deeply RNA can influence sequence representations.

**Recommendation:** Start with approach (b). It is fastest to iterate on and the Flamingo paradigm is well-validated. If frozen cross-attention proves insufficient, move to (a) or continued pretraining at moderate scale.

---

## 5. Adapter / Prefix Approaches — Lightweight Conditioning

### Adapter layers (Houlsby et al. 2019)
- Insert small bottleneck MLP modules after each Transformer layer.
- The RNA signal is encoded separately and added to the adapter input.
- Very parameter-efficient (~0.5–2% of base model params).

### Prefix tuning / prompt tuning (Li & Liang 2021)
- Prepend a learned "prefix" of virtual tokens to the input sequence.
- These prefix tokens can be conditioned on RNA embeddings (soft prompts computed from RNA data).
- Almost no parameter overhead.
- Limitation: the prefix can only influence the sequence through attention from sequence to prefix — weaker than cross-attention injected at every layer.

**Genomic application:** Encode the RNA tracks for a given locus into a small set of "context tokens" (e.g., 64 tokens representing the expression state of this cell). Prepend these to the DNA sequence tokens before passing through the frozen DNA LM. Analogous to a "system prompt" for sequence context.

### LoRA (Hu et al. 2021)
- Low-rank decomposition of weight matrix updates: ΔW = AB.
- Allows updating Evo2 (7B params) weights at very low cost.
- Can be combined with cross-attention adapters.

### FiLM conditioning (Perez et al. 2018)
- Encode **which cell type / condition** as a discrete embedding, then inject via Feature-wise Linear Modulation at each layer:
  `h_layer = LayerNorm(h) * gamma(cell_embedding) + beta(cell_embedding)`
- Exactly what Enformer/Borzoi do for multi-track prediction.
- Extremely cheap. Fails for continuous variation across genotypes or novel cell types not seen at training.

---

## 6. Graph-Based Approaches — Gene Regulatory Networks as Structure

### The idea

Instead of treating the genome as a linear sequence, represent it as a graph where:
- **Nodes** = genes, enhancers, regulatory elements
- **Edges** = regulatory interactions (TF binding, chromatin loops, eQTLs, co-expression)

RNA expression values are node features. DNA sequence features (from the LM) are also node features. A GNN propagates information along regulatory edges.

### Relevant models

**SCENIC / GRNBoost2 (Aibar et al. 2017):** Infers gene regulatory networks from scRNA-seq. The resulting GRN can define the graph structure for a GNN.

**GEARS (Roohani et al. 2023, Nature Biotechnology):** GNN-based perturbation effect prediction using gene co-expression graph. Directly relevant for learning gene covariation.

**Chromatin loop / Hi-C graph:** Model the genome as a graph where edges come from Hi-C contact frequency. Node features: DNA sequence (from frozen LM embedding), RNA coverage. Message passing propagates RNA expression signals along chromatin contacts.

**CellOracle (Kamimoto et al. 2023, Nature):** Transcription factor network model using scRNA + sequence motifs.

### Tradeoffs

| Advantage | Limitation |
|---|---|
| Naturally handles long-range regulatory interactions | Graph structure must be pre-defined or learned separately |
| Can incorporate prior knowledge (known regulatory networks) | GRNs are noisy and cell-type-specific |
| Scalable to genome-wide analysis | Loses base-pair resolution; works at gene/element level |
| Strong inductive bias for co-regulation | Hard to train end-to-end with a DNA LM |

### Hybrid: GNN + DNA LM

1. DNA LM embeds each regulatory element (gene body, enhancer) → node embedding.
2. RNA coverage for each element → additional node feature.
3. Hi-C or co-expression edges define the graph.
4. GNN aggregates neighboring regulatory signals.

---

## 7. What Architecture Choices Matter Most for Learning Gene Covariation

### 7.1 Column attention across samples (most critical)

The MSA Transformer's core finding: to learn co-evolutionary constraints, you need attention **across sequences at the same position** (column attention). The genomic analogue is attention across cells/genotypes at the same genomic coordinate. Without this, the model cannot learn that "gene A is high when gene B is high" because it never jointly attends to both.

Column attention is cheap when M (number of cells) is small. Even with M = 50–200 cells, column attention is efficient (M² is small).

### 7.2 Multi-cell batching with explicit genotype labels

Covariation is a property of **sets of cells**, not individual cells. Training batches should be "groups" of cells from the same locus, different genotypes. The loss should reward predicting expression of one cell given the others — like masked language modeling over the cell axis (analogous to how the MSA Transformer masks random rows and predicts from remaining ones).

### 7.3 Positional encoding at genomic scale

Gene covariation is partly positional: co-regulated genes tend to be in the same TAD, share enhancers, or are physically co-localized.
- **RoPE / ALiBi:** Distance between tokens decays attention weight. RoPE is used in Evo2 and generalizes to long contexts.
- **Enformer-style binning:** Group positions into 128 bp bins; attention at bin level captures TAD-scale interactions.

### 7.4 Contrastive / self-supervised objective over cells

A powerful training signal: given cells from two different genotypes at the same locus, the model should learn representations that differ in ways predictable from the DNA sequence differences. Analogues: DNABERT-S (contrastive learning on DNA), scVI (VAE with explicit genotype latent).

### 7.5 Sparse attention for long genomic contexts

At base-pair resolution, L = 200k for a gene body. Options:
- **Sliding window attention (Longformer):** Each position attends to a local window.
- **Hierarchical processing:** CNN/pooling to 128 bp → full attention at 128 bp resolution. Most practical choice.
- **Linear attention (Hyena):** O(L) complexity. Already in Evo2. Natural choice for the RNA encoder.

### 7.6 Explicit modeling of zero-inflation

Single-cell RNA-seq has high dropout:
- **Zero-inflated negative binomial (ZINB) output distribution** (as used in scVI, scANVI).
- **Masking:** Mask zero positions during column attention.

---

## 8. Less Obvious Approaches Worth Considering

### 8.1 Genotype-conditioned expression prediction as the primary objective

Rather than passing RNA as input, train Evo2/NTv3 to **predict** per-genotype RNA coverage from the individual's sequence. Use Smart-seq from many genotypes as training data. This is the "personalized Enformer" direction (Huang et al. 2023, bioRxiv) and may be the most direct path to your goal — it uses your Smart-seq data as a training signal, leverages pretrained DNA LMs, and explicitly forces learning of genotype → expression covariation.

### 8.2 Retrieval-augmented DNA models

Given a query DNA sequence, retrieve the most similar genotypes from a database of (genotype, RNA track) pairs and condition on retrieved tracks via cross-attention. Analogues: RETRO (Borgeaud et al. 2022), Atlas (Izacard et al. 2023). Avoids needing all genotypes simultaneously at training time.

### 8.3 Causal perturbation data as covariation supervision

Perturbation of gene A (CRISPRi/a) changing expression of gene B is **causal** covariation. Relevant datasets: **Genome-wide Perturb-seq (Replogle et al. 2022, Cell)**, **sci-Plex**. Adding a perturbation-prediction objective forces learning of regulatory logic rather than just correlation.

### 8.4 Splicing-aware token design

Smart-seq full-coverage data contains genuine splicing information. Do not pool to gene-level counts before encoding. Keep coverage as a function of position. 128 bp binning preserves exon-level features. SpliceAI (Jaganathan et al. 2019) shows dilated CNNs with ~10k bp receptive fields capture exon-intron boundary signals.

### 8.5 RNA velocity as a temporal signal

RNA velocity (La Manno et al. 2018; Bergen et al. 2020, scVelo) estimates the time-derivative of gene expression from spliced/unspliced ratios — directly available from Smart-seq full-coverage data (intronic reads = unspliced). Use velocity as an additional channel alongside steady-state coverage to provide a directional/temporal signal.

### 8.6 Multi-species sequence alignments as a free MSA

Use **multi-species sequence alignments** (100-way UCSC vertebrate alignment, Zoonomia 240-mammal alignment) as the MSA input before investing in RNA conditioning infrastructure. Evolutionary covariation is strong signal for regulatory function — AlphaFold2 used this exact strategy. This is a compelling baseline.

---

## 9. Synthesized Recommendations

### Tier 1: Highest impact, most tractable

1. **Genotype-conditioned expression prediction** (supervised, Enformer-style but with personalized sequences): Train Evo2/NTv3 to predict per-cell-type RNA-seq coverage from the individual's sequence. The covariation signal comes from the training distribution. Most direct path.

2. **MSA-style axial attention with 128 bp tokenization**: Joint model with row attention (within each cell's RNA track) and column attention (across cells at the same genomic position). Use convolutional front-end to reduce to 128 bp tokens. Train with masked prediction over the cell axis. Directly operationalizes the researcher's intuition.

### Tier 2: Strong alternatives worth prototyping

3. **Flamingo-style frozen DNA LM + RNA cross-attention**: Freeze Evo2, train a Hyena-based RNA encoder, insert cross-attention layers. Fastest path to a working prototype.

4. **Multi-species MSA baseline**: Use the Zoonomia 240-mammal alignment as an MSA input before building RNA infrastructure. Validates the axial attention setup at lower data acquisition cost.

### Tier 3: Longer-term / higher-risk

5. **Graph + DNA LM hybrid**: Hi-C-based regulatory graph with DNA LM embeddings as node features.

6. **Perturbation-supervised joint training**: Add Perturb-seq data as a causal covariation signal.

---

## 10. Key Papers

| Paper | Relevance |
|---|---|
| Rao et al. 2021, "MSA Transformer" (NeurIPS) | Core MSA axial attention architecture |
| Jumper et al. 2021, "AlphaFold2" (Nature) | Full MSA stack + pair representation + triangular updates |
| Avsec et al. 2021, "Enformer" (Nature Methods) | Multi-track expression prediction from DNA; 128 bp tokenization |
| Linder et al. 2023, "Borzoi" (Nature Genetics) | RNA-seq track prediction at bp resolution; U-Net decoder |
| Alayrac et al. 2022, "Flamingo" (NeurIPS) | Cross-modal attention, conditioning frozen LM on new modality |
| Theodoris et al. 2023, "Geneformer" (Nature) | Large-scale scRNA pretraining; cell state representations |
| Wang et al. 2023, "scGPT" (Nature Methods) | GPT-style pretraining on scRNA; perturbation prediction |
| Roohani et al. 2023, "GEARS" (Nature Biotechnology) | GNN for perturbation-driven gene covariation |
| Nguyen et al. 2024, "Evo2" (Arc Institute) | Base DNA LM; StripedHyena architecture |
| Bergen et al. 2020, "scVelo" (Nature Biotechnology) | RNA velocity from spliced/unspliced; temporal signal |
| Jaegle et al. 2021, "Perceiver IO" (ICML) | Cross-modal compression via cross-attention |
| Zhou 2022, "Orca" (Nature Methods) | Hi-C contact prediction from sequence |
| Replogle et al. 2022, Genome-wide Perturb-seq (Cell) | Large-scale causal covariation dataset |
| Jaganathan et al. 2019, "SpliceAI" (Cell) | Splicing prediction from sequence; dilated CNN architecture |
| Lopez et al. 2018, "scVI" (Nature Methods) | Probabilistic model for scRNA; ZINB, batch correction |
| Perez et al. 2018, "FiLM" (AAAI) | Feature-wise linear modulation for multi-task conditioning |

---

*Notes: Knowledge cutoff August 2025. Check bioRxiv for 2025 preprints on "multi-modal genomics language model", "single-cell DNA language model", and "genotype-expression model" for the latest developments.*
