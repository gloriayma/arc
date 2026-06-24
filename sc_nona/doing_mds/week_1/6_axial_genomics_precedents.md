# Axial / Stack-Attention Precedents in Genomics

Scope: precedents for an Evoformer-style architecture over a *variable-depth* stack of *aligned, base-pair-resolution* genomic tracks (e.g. smart-seq scRNA coverage across cells) plus DNA, with attention along both the genome axis and the stack axis.

Depth markers used below: **[read closely from search]** = read paper/abstract carefully via search summary, **[skimmed]** = read search summary only, **[name only]** = aware of, did not verify details.

## TL;DR

- **The exact architecture proposed — Evoformer-style row+column axial attention over a stack of aligned base-pair coverage tracks where stack depth varies per training example and tracks come from scRNA-seq of different cells — does not appear to exist in the published genomics literature as of mid-2026.** It is a real gap.
- The **closest single precedent is GPN-MSA / GPN-Star** (Benegas et al., Songlab). They run a transformer over a stack of *aligned DNA sequences across vertebrate species* with masked-LM. That validates "aligned stack + axial-ish attention works on genomes," but the stack is DNA, not RNA expression, and GPN-MSA appears to *collapse the species dimension into a per-column feature vector* rather than doing true MSA-Transformer-style axial attention.
- The **closest single-cell precedent is scooby** (Gagneur lab, Nat. Methods 2025): per-cell scRNA + scATAC coverage from sequence using a Borzoi backbone + a *cell-conditioned decoder*. scooby does NOT do attention across cells — the cell axis is a *condition*, not a *set to attend over*.
- **Decima** (Genentech / gReLU) is similar to scooby but at pseudobulk resolution: sequence -> per-pseudobulk-track expression vector. Again, no cross-sample attention; the sample axis is the output channel.
- **No work treats a variable-depth set of base-pair coverage tracks as an MSA-style input to be attended over, on either axis.** The proposed architecture sits in a real white space at the intersection of (a) AF-style axial attention on biological stacks, (b) sequence-to-track regulatory models, and (c) sc / full-length RNA-seq data.

## 1. Closest precedents in genomics

### 1.1 GPN-MSA — alignment-based DNA LM over a stack of species [read closely from search]
- Benegas, Albors, Aw, Ye, Song. *A DNA language model based on multispecies alignment predicts the effects of genome-wide variants.* Nat. Biotechnol. 2025 (bioRxiv 2023-10-10). https://www.biorxiv.org/content/10.1101/2023.10.10.561776 ; https://www.nature.com/articles/s41587-024-02511-w
- Input: a 128-bp window of a 100-vertebrate whole-genome MSA. Columns = positions on the human genome; rows = aligned species.
- The species dimension is encoded into a per-column feature vector and a standard 1D transformer runs across columns. The "rows" are summarized, not attended over. (Worth verifying against the repo's existing `gpn_arch.html`.)
- Trained with masked-LM on the human row. SOTA on ClinVar / COSMIC / OMIM / gnomAD; beats CADD / phyloP / Enformer / SpliceAI. 3.5 hours on 4× A100.
- **Validates vertical stacking + masked-LM as a pretraining signal.** Adding true row+column axial attention here is essentially open work.

### 1.2 GPN-Star — phylogeny-aware MSA model [skimmed]
- Benegas et al. *Predicting functional constraints across evolutionary timescales with phylogeny-informed genomic language models.* bioRxiv 2025-09-21. https://www.biorxiv.org/content/10.1101/2025.09.21.677619
- Successor to GPN-MSA. Adds species-tree structural prior on the row axis. Vertebrate / mammal / primate timescales; also mouse, chicken, fly, worm, Arabidopsis.
- Closest precedent for *non-trivial structure on the stack axis*. Analog for scRNA: cells aren't exchangeable — they have a cell-type / lineage tree.

### 1.3 PhyloGPN — phylogeny in the loss, not the architecture [skimmed]
- arXiv:2503.03773 / PMC11908359, 2025.
- CNN maps a 481-bp single human DNA window to F81 substitution-model parameters; Felsenstein's pruning algorithm computes alignment likelihood on the species tree; NLL loss.
- Phylogeny enters only the loss; inference is single-sequence (alignment-free).
- Contrasting design point: embed evolutionary info into the *training signal* rather than the input. The scRNA analog is essentially what scooby/Decima do.

### 1.4 scooby — per-cell sc-omics coverage from sequence [read closely from search]
- Yang, Hingerl, Eraslan, ... Gagneur. *scooby: modeling multimodal genomic profiles from DNA sequence at single-cell resolution.* Nat. Methods 2025 (bioRxiv 2024-09-19). https://www.biorxiv.org/content/10.1101/2024.09.19.613754 ; https://github.com/gagneurlab/scooby
- Borzoi backbone (parameter-efficient fine-tune) at 32-bp resolution; small *cell-state-specific decoder* takes (sequence embedding, cell embedding) -> scRNA-seq coverage and scATAC insertion profile.
- **Cell axis = condition, not set.** No attention across cells. For N cells you run the decoder N times. Backbone is cached across cells.
- Direct baseline for the proposed approach. Weakness: no covariation modeling inside the model.
- The repo already has `doing_mds/3.31_scooby_borrowables.md`, so this is on the user's radar.

### 1.5 Decima — pseudobulk sequence-to-expression [skimmed]
- Lal, Wenkel, ... Avsec. *Decoding sequence determinants of gene expression in diverse cellular and disease states.* bioRxiv 2024-10-09. https://www.biorxiv.org/content/10.1101/2024.10.09.617507
- 22M cells aggregated into 8,856 pseudobulks × 18,457 genes. Borzoi-like backbone; output is a vector over pseudobulks.
- No attention across samples. Reinforces the prevailing pattern: sample axis is output or condition, never input-to-attend-over.

### 1.6 Enformer / Borzoi / EpiBERT / gReLU — tracks as output [read closely from search]
- Enformer (Avsec 2021); Borzoi (Linder, Srivastava, Kelley, Nat. Genet. 2024, https://www.nature.com/articles/s41588-024-02053-6); EpiBERT (Javed et al., Cell Genomics 2025, https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00018-7); gReLU framework (Lal et al., Nat. Methods 2025, https://www.nature.com/articles/s41592-025-02868-z).
- DNA is the only sequential input. Sample / cell identity enters via output channel index or cell-embedding condition. None do attention over a stack of per-sample tracks as input.

### 1.7 EpiFoundation / ChromFound / GET / EpiAttend [skimmed]
- EpiFoundation: bioRxiv 2025-02-05, https://www.biorxiv.org/content/10.1101/2025.02.05.636688 . scATAC foundation model via peak-to-gene alignment.
- All operate at peak / gene granularity, not base-pair. Same family as scGPT / Geneformer with chromatin tokens.

### 1.8 Cisformer — cross-modality (not cross-cell) attention [skimmed]
- Genome Biology 2025. https://link.springer.com/article/10.1186/s13059-025-03823-z . Cross-attention between RNA and ATAC modalities of the same cell. Axis is modality, not cell.

## 2. Adjacent work informing the design

### 2.1 NTv3 — single-bp, 1 Mb context, species-conditioned [skimmed]
- InstaDeep, Dec 2025. https://instadeep.com/2026/02/modelling-the-genome-with-ntv3/ ; https://instadeep.com/wp-content/uploads/2025/12/NT_v3.pdf
- U-Net (conv-down + transformer core + conv-up), single-bp tokenization, 1 Mb context.
- **Multi-species via adaptive layer norm on species token**, not stacking. Direct alternative design point to your proposed stack-attention.
- Post-trained on ~16k functional tracks from 24 species.

### 2.2 Borzoi-derived personal-genome models [skimmed]
- Borzoi fine-tuned on 54,219 UKB individuals × 2,923 plasma proteins (bioRxiv 2025-09-26, https://www.biorxiv.org/content/10.1101/2025.09.26.678908). Per-gene per-individual models; +86% VEP improvement.
- "Personal genome via fine-tune," not "personal genomes as input set."

### 2.3 ChatNT / DNA-Perceiver — Perceiver-style for DNA [skimmed]
- ChatNT (Nat. Mach. Intell. 2025, https://www.nature.com/articles/s42256-025-01047-1): Flamingo Perceiver Resampler with gated cross-attention bridges NT embeddings into an LLM.
- DNAPerceiver / ProteinPerceiver (Comp Methods Programs Biomed 2023, https://pubmed.ncbi.nlm.nih.gov/37004267/).
- The variable-size set being compressed is the sequence, not a stack of per-sample tracks. Still the right toolkit for variable-depth handling.

### 2.4 Boltz-1 / Chai-1 / Protenix — efficient MSA depth [skimmed]
- Boltz-1 (Wohlwend et al., 2024, https://pmc.ncbi.nlm.nih.gov/articles/PMC11601547/): Dense MSA Pairing Algorithm; chunking in transition / pair-weighted-average / outer-product layers.
- Chai-1: runs with or without MSA at near-equal performance.
- Protenix-Mini+ (arXiv:2510.12842): scalable pairformer.
- Variable-depth MSAs have been made practical via chunking + sparse pairing. Directly transferable engineering ideas.

### 2.5 scGAA — axial attention in scRNA, wrong axes [skimmed]
- Sci. Rep. 2024. https://www.nature.com/articles/s41598-024-73356-1
- Axial attention over cells × genes count matrix. Not bp × cells.

### 2.6 Selective SSMs for RNA-seq coverage [skimmed]
- bioRxiv 2025-02-13, https://www.biorxiv.org/content/10.1101/2025.02.13.638190
- Mamba-style SSMs beat transformers at Borzoi-style coverage prediction on the genome axis. Suggests: genome axis may not need attention; use attention only on cell axis.

### 2.7 Pangenome-aware methods [skimmed]
- Pangenome-aware DeepVariant (bioRxiv 2025-06-05) etc. Not transformer-over-stack; no precedent for axial attention over individuals.

## 3. Honest gap assessment

**What does not exist (best-effort search):**
1. Evoformer-style axial attention over a variable-depth stack of base-pair-resolution coverage tracks from sc / bulk RNA-seq.
2. Treating cells as a permutation-equivariant set attended over jointly with the genome axis.
3. Cross-gene expression covariation learned via attention on aligned bp-resolution coverage tracks (closest = scGPT/scFoundation on count matrices).
4. True row+column MSA-Transformer-style axial attention adapted to genome alignments — surprising gap given how natural the port from MSA Transformer would be.

**What exists, should not be re-invented:**
1. Per-cell decoders on a shared sequence backbone (scooby, Decima) — your baseline.
2. Stacked-species masked-LM (GPN-MSA, GPN-Star) — useful as initialization or DNA-only control.
3. Variable-depth-stack engineering (Boltz-1 chunking, Perceiver-resampler compression).

**What differentiates the proposed approach:**
- Input axis 1 = genome position. Input axis 2 = variable set of cells with aligned bp coverage tracks. Rows carry *functional* (transcription) signal, not just *evolutionary* DNA.
- Axial attention enables learning gene-gene covariation across cells AND long-range promoter-enhancer pairing along the genome.
- Permutation equivariance on cell axis — no fixed cell ontology, any depth.
- Variation across genotypes on cell axis forces DNA -> RNA causal modeling.

Closest neighbors: GPN-MSA (right axial structure, wrong row signal) and scooby (right row signal, wrong axial structure). Filling the diagonal between them is genuine white space.

## 4. Open questions

1. **Row encoder for a coverage track.** Bp-scalar coverage tokenization (binning? per-bp linear? small conv?) governs whether bp resolution survives stack-attention.
2. **DNA-anchor vs MSA-rep pairing.** Is DNA the "single representation" (Evoformer terminology) and coverage tracks the "MSA representation," or is DNA an extra row? AF3-style single↔MSA interaction matters.
3. **Cell-axis structure.** GPN-Star says rows aren't exchangeable. Inject cell-type / donor / tissue via biased attention or tree-conditioned position rather than pure set.
4. **Why hasn't this been done?** (a) Full-length scRNA with population-scale genotype variation is scarce — Smart-seq2/3 + diverse genotypes is the data bottleneck. (b) Memory: bp × cells tensors blow up; Boltz-1 chunking needed. (c) Sociological: regulatory community = tracks-as-output (Enformer line); single-cell community = counts-as-tokens (scGPT line). The cross-pollination required is non-obvious.
5. **Splicing from bp coverage?** Junction-spanning reads create discontinuities encoding splice events implicitly. But may need explicit junction tracks (Borzoi does) — worth empirically checking.
6. **Column attention at 100 kb+.** At 1 Mb × single-bp × hundreds of cells, axial attention is feasible only with strong sparsification on the genome axis (windowed / Mamba-on-genome, attention-on-cells). Perceiver-resampler variant — compress cell stack to small latent set then cross-attend from genome — may be a more practical first version than full Evoformer.
7. **Initialization.** Backbone from NTv3 or Borzoi; cell-axis attention from scratch. Most practical training plan.

## 5. Selected references

- GPN-MSA, Nat. Biotechnol. 2025 — https://www.nature.com/articles/s41587-024-02511-w **[read closely from search]**
- GPN-Star, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.09.21.677619 **[skimmed]**
- PhyloGPN, arXiv:2503.03773 **[skimmed]**
- scooby, Nat. Methods 2025 — https://www.biorxiv.org/content/10.1101/2024.09.19.613754 **[read closely from search]**
- Decima, bioRxiv 2024 — https://www.biorxiv.org/content/10.1101/2024.10.09.617507 **[skimmed]**
- Borzoi, Nat. Genet. 2024 — https://www.nature.com/articles/s41588-024-02053-6 **[name only]**
- EpiBERT, Cell Genomics 2025 — https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00018-7 **[skimmed]**
- EpiFoundation, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.02.05.636688 **[name only]**
- Cisformer, Genome Biol. 2025 — https://link.springer.com/article/10.1186/s13059-025-03823-z **[name only]**
- NTv3, InstaDeep Dec 2025 — https://instadeep.com/wp-content/uploads/2025/12/NT_v3.pdf **[skimmed]**
- Boltz-1, 2024 — https://pmc.ncbi.nlm.nih.gov/articles/PMC11601547/ **[skimmed]**
- gReLU, Nat. Methods 2025 — https://www.nature.com/articles/s41592-025-02868-z **[name only]**
- ChatNT, Nat. Mach. Intell. 2025 — https://www.nature.com/articles/s42256-025-01047-1 **[skimmed]**
- Selective SSMs for RNA-seq coverage, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.02.13.638190 **[name only]**
- scGAA, Sci. Rep. 2024 — https://www.nature.com/articles/s41598-024-73356-1 **[name only]**
