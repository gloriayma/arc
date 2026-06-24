# CURRENT_LIMS.md — Limitations and Failure Modes of DNA/Genomic Language Models

*Compiled from training knowledge (cutoff August 2025). Performance numbers and benchmark details are drawn from training memory; treat specific figures as approximate and verify against primary literature before citing.*

---

## 1. Benchmarks That Exist — and Where Models Underperform

### 1.1 Nucleotide Transformer Benchmarks (NT-v1/v2/v3)

The Nucleotide Transformer papers (Dalla-Favera et al., INRIA/EMBL) introduced a suite of 18 downstream tasks used to compare models. These include:

- **Histone modification prediction** (H3, H3K4me1, H3K4me2, H3K4me3, H3K9ac, H3K14ac, H3K36me3, H3K79me2, H4, H4ac) — short 200 bp windows, essentially a sequence→binary classification. Models generally perform well here (MCC > 0.7 for most marks in NT-v2), because the signal is locally encoded.
- **Promoter prediction** (TATA vs. non-TATA, all promoters) — again short-range. Models do well.
- **Splice site prediction** (donor/acceptor) — moderate performance; models can learn the canonical GT-AG rule but struggle on weak or tissue-specific sites.
- **Enhancer prediction** — binary classification on ENCODE-derived sequences. Performance is decent for strong enhancers but drops sharply for cell-type-specific or developmental enhancers.

**Where models underperform on NT benchmarks:** The benchmarks are mostly binary classification on fixed short windows. They do not test *functional specificity* — a model can predict "this is an enhancer" without understanding which gene it regulates or in which cell type.

---

### 1.2 Genomic Benchmarks (Grešová et al., 2023)

A standardized collection of 8 classification tasks including mouse non-coding RNA, human nontata promoters, human enhancers (Cohn lab), regulatory sequences, and splice sites. Used to compare DNABERT, HyenaDNA, Nucleotide Transformer, etc.

**Key finding from the benchmark paper:** No model dominates across all tasks. HyenaDNA's long-context advantage does not always translate to better performance on tasks where the relevant signal is local. DNABERT (k-mer tokenized, BERT-style) remains competitive on short-range tasks.

**Weakness exposed:** The benchmark does not include any tasks that require understanding *relationships between genomic elements* — it is purely element-level classification.

---

### 1.3 BEND Benchmark (Marin et al., 2023 — Berman lab)

BEND (Benchmark for DNA Language Models) is more comprehensive and includes:

- Gene finding (multi-class: coding, non-coding, intergenic)
- Chromatin accessibility (ATAC-seq peak prediction)
- Histone modification
- Transcription factor binding
- Variant effect prediction (eQTL-linked SNPs)
- Enhancer-gene association (a coarser version)

**Reported gaps in BEND:** Variant effect prediction tasks are where most models show the largest performance gaps relative to task-specific supervised models. Models trained purely on sequence struggle to predict the *direction* and *magnitude* of expression changes from SNPs, even common ones. The eQTL task in BEND tests whether a model can rank SNPs by predicted effect size — DNA LMs perform only modestly above random on this.

---

### 1.4 Enformer and Successor Models

Enformer (Avsec et al., Nature Methods 2021, DeepMind) is not a language model in the NLP sense — it is a convolutional + transformer model trained with supervised labels (CAGE, DNase, ChIP-seq, RNA-seq tracks from hundreds of cell types). It accepts 196 kb of sequence and predicts ~5,000 genomic tracks.

**What Enformer is actually good at:**
- Predicting bulk-level gene expression from sequence in well-studied tissues (Pearson r ≈ 0.85–0.90 on held-out genes for common cell types)
- Predicting effects of non-coding variants on gene expression (r ≈ 0.6–0.7 on GTEx eQTLs, depending on tissue)

**Where Enformer fails:**
- **Rare cell types and developmental contexts:** Trained on bulk ENCODE/FANTOM data, it has no representation of cell-type-specific regulation that emerges only in single-cell data.
- **Long-range beyond 196 kb:** Enformer's receptive field is ~196 kb. Enhancer-gene pairs beyond ~100 kb from the TSS are systematically underweighted. Hi-C and CTCF loop data are not in the training signal.
- **Variant effect saturation:** Enformer-based in silico saturation mutagenesis (ISM) loses specificity in repeat-rich regions and at positions far from annotated regulatory elements.
- **Cross-species generalization:** Performance on non-human species is substantially worse unless fine-tuned on species-specific data.

The Borzoi model (Linder et al., 2023) extended Enformer to ~524 kb and added RNA-seq tracks (not just CAGE), yielding improvements on splicing and isoform-level predictions. Still does not incorporate single-cell resolution.

---

### 1.5 Evo2 (Arc Institute, 2024–2025)

Evo2 is a 7B (and 40B) parameter model trained on ~9 trillion nucleotides across diverse organisms using the StripedHyena architecture (a hybrid of Mamba-style SSMs and attention). Its context window is 1 million base pairs.

**What Evo2 claims to do well:**
- Zero-shot fitness prediction for proteins and non-coding regulatory elements
- Variant effect prediction across species
- Generation of plausible synthetic genomic sequences

**Reported limitations:**
- Evo2's evaluations are mostly on *functional variant effects in model organisms* (yeast, bacteria) or on protein-coding sequences where ground truth fitness is available from DMS experiments. Human non-coding regulatory evaluation is limited.
- The model is evaluated on splice site prediction and gains from long context — but this is a local signal task, so gains are modest.
- **Promoter-enhancer pairing is not a direct evaluation task in Evo2's paper.** The model can presumably learn co-occurrence statistics in genomic sequence, but there is no benchmark showing it can identify which enhancer drives which gene.
- Long-context utilization: the 1M bp window is architecturally available, but ablations in the Evo2 paper show that for most tasks, most of the gain comes from ~10–50 kb context. Utilization of the full million-base context is not well demonstrated.

---

### 1.6 Caduceus (Schiff et al., 2024)

Caduceus uses a bidirectional Mamba (RC-equivariant by construction) for sequence modeling. Evaluated on Genomic Benchmarks and Long Range Arena (genomics variant).

**Key finding:** Caduceus outperforms HyenaDNA on most Genomic Benchmarks tasks, particularly those requiring long-range context (enhancer-promoter orientation tasks). However, the long-range tasks in these benchmarks are still relatively short (sequences up to ~1–4 kb). True multi-kilobase enhancer-promoter pairing is not tested.

---

### 1.7 Variant Effect Prediction — ClinVar / gnomAD Tasks

Multiple groups have benchmarked DNA LMs on ClinVar pathogenicity classification and gnomAD allele frequency prediction (as a proxy for constraint).

- **CADD, SIFT, PolyPhen** remain strong baselines for coding variants; DNA LMs close the gap on non-coding variants but typically underperform task-specific tools like DeepSEA, Sei, or Enformer for regulatory variant effects.
- **Reported numbers (approximate, from training memory):** NT-v2 on ClinVar classification achieves ~70–75% accuracy; task-specific models reach ~80–85%.

---

## 2. Tasks Models Should Be Evaluated On But Aren't

### 2.1 Enhancer-Gene Pairing

**The core gap the researcher identified.** No standard benchmark directly tests whether a model can predict, given an enhancer sequence and a gene's promoter sequence, whether they interact. Existing enhancer benchmarks test *is this sequence an enhancer* — not *which gene does this enhancer regulate*.

Partial proxies exist (Hi-C-based loop prediction, Activity-by-Contact scores) but are not in standard DNA LM evaluation suites. The ABC model (Fulco et al., 2019) provides ground-truth enhancer-gene links from CRISPRi screens, but this has not been formalized as a language model benchmark.

**Why this is hard without RNA supervision:** The model needs to learn that enhancer A and gene B co-vary in expression across cell types, which is a signal that only emerges in multi-condition RNA data. Sequence alone encodes the *potential* for interaction, not its *realization* in a specific cell state.

### 2.2 Cell-Type-Specific Regulatory Grammar

Current benchmarks evaluate on bulk-averaged ENCODE data. There is no standard benchmark asking: "Given this regulatory sequence, predict its activity specifically in hepatocytes vs. cardiomyocytes." scRNA-seq + scATAC-seq multiome data provides exactly this supervision — but it has not been systematically turned into DNA LM benchmarks.

### 2.3 Isoform-Level Prediction in Specific Cell Types

Models like Borzoi predict aggregate isoform ratios, but no benchmark tests prediction of *cell-type-specific splicing* from sequence. Long-read single-cell data (PacBio MAS-seq, Oxford Nanopore) is starting to provide ground truth for this.

### 2.4 cis-Regulatory Module Composition

Can a model recognize that two transcription factor binding sites in a specific spatial arrangement create a synergistic regulatory module? MPRA (massively parallel reporter assay) datasets encode this, but they are not standard DNA LM benchmarks.

### 2.5 3D Genome Organization

No standard benchmark tests whether DNA LMs can predict TAD boundaries, compartment identity (A/B), or loop anchor positions from sequence. Orca (Zhou et al., 2022) does this with a supervised approach, but it is not used to evaluate foundation models. TAD boundaries have strong sequence signatures (CTCF motifs, cohesin binding sites) — there is real opportunity here that is not exploited.

### 2.6 Temporal / Developmental Trajectories

Can a model predict how regulatory accessibility changes during differentiation? This requires connecting sequence features to dynamic RNA/ATAC profiles, a supervision signal only available from trajectory-resolved single-cell data. No such benchmark exists for DNA LMs.

### 2.7 Cross-Cell-Type Transfer of Variant Effects

A given non-coding variant may have different effects in different cell types (cell-type-specific eQTLs from GTEx or OneK1K). Current variant effect benchmarks aggregate across tissues. A benchmark that tests cell-type-stratified variant effects would be more informative and does not currently exist in DNA LM evaluation.

---

## 3. Long-Context Limitations

### 3.1 Architecture vs. Utilization

Most DNA LMs advertise large context windows but do not demonstrate proportional gains from using them.

| Model | Context window | Demonstrated useful range (approx.) |
|---|---|---|
| DNABERT-2 | 512 bp (BPE tokenized) | 512 bp |
| HyenaDNA | Up to 1M bp | ~1–10 kb for most tasks |
| Nucleotide Transformer v2 | 6 kb | 6 kb |
| Caduceus | ~4 kb (benchmarked) | ~1–4 kb |
| Enformer | 196 kb (fixed) | ~50–100 kb |
| Borzoi | ~524 kb | ~200 kb |
| Evo2 | 1M bp | ~10–50 kb for demonstrated tasks |

The gap between architectural context and *demonstrated useful context* is large. This is partly a training data issue: the pretraining objective (next-token prediction or MLM on genomic sequence) does not force the model to use long-range information unless long-range dependencies are necessary to predict the masked/next token. For most of the genome, local context is sufficient to predict the next nucleotide.

### 3.2 Positional Encoding Degradation

Transformer-based DNA LMs use learned positional embeddings (DNABERT, NT) or ALiBi/RoPE variants. Extrapolation beyond the training context length is unreliable. SSM-based models (HyenaDNA, Mamba-based Caduceus, Evo2's StripedHyena) have better theoretical scaling, but in practice the long-range signal is not learned unless the pretraining tasks reward it.

### 3.3 The Enhancer-Promoter Distance Problem

The median distance from an enhancer to its target gene TSS in humans is ~50–200 kb (based on HiChIP and CRISPRi-distance studies). The 95th percentile extends beyond 500 kb. This means:

- Enformer (196 kb) misses a significant fraction of real enhancer-gene pairs
- Models with smaller windows (NT, DNABERT, Caduceus at standard benchmarked lengths) miss the vast majority
- Even Evo2's 1M bp window *architecturally* covers most pairs, but whether it *learns* to use this information is undemonstrated

### 3.4 Repeat Masking and Data Quality

A large fraction of the human genome is hard-masked repeats in training data. Models trained on the reference genome with N-masking for repeats have explicit blind spots in ~50% of the genome. This affects heterochromatic gene regulation, LINE/SINE-derived regulatory elements (which are real and important), and telomeric/centromeric sequences.

---

## 4. Gene-Level Understanding

### 4.1 Do Models Know What a Gene Is?

**Short answer: weakly, and in a shallow way.**

DNA LMs can learn to predict splice sites and coding sequence boundaries from sequence context, which implicitly encodes gene structure knowledge. Probing experiments have not been extensively published for DNA LMs, but available evidence suggests:

- Models trained on raw sequence develop internal representations that correlate with exon/intron identity and CpG island content
- Gene body methylation patterns and codon usage biases are likely captured implicitly
- But there is no evidence that models have learned a *symbolic* gene-level representation — they operate at the nucleotide/token level

### 4.2 Gene Identity vs. Gene Function

A DNA LM, given the sequence of two paralogs (e.g., BRCA1 and a pseudogene), has no mechanism to distinguish them by *function* unless functional differences are encoded in sequence divergence. The model cannot reason about protein domains, PPI networks, or pathway membership.

### 4.3 Gene-Gene Relationships

DNA LMs have essentially no capacity to model *trans*-regulatory relationships, because these relationships are not encoded in the linear sequence of either gene. What is missed:

- Transcription factor → target gene regulatory cascades
- microRNA → mRNA suppression (requires knowing both sequences AND their expression context)
- Polycomb domain silencing (requires 3D genome structure knowledge)
- Synthetic lethality relationships

### 4.4 What RNA Supervision Could Fix

If a DNA LM is jointly trained or fine-tuned with scRNA-seq data providing cell-type-specific expression profiles, it gains access to:

1. **Co-expression structure** — genes co-regulated across cell types likely share regulatory grammar
2. **Promoter-enhancer pairing signals** — if enhancer accessibility (scATAC-seq) co-varies with gene expression (scRNA-seq) across cell types, this pair-wise co-variation provides indirect supervision for enhancer-gene pairing
3. **Dosage sensitivity** — genes that vary a lot across cell types vs. housekeeping genes have different regulatory architectures
4. **Developmental trajectories** — expression dynamics along pseudotime constrain which regulatory elements are causal

---

## 5. Surprising Findings and Outside-the-Framing Observations

### 5.1 Models Are Better at Coding Sequence Than Often Assumed

DNA LMs implicitly learn codon usage, reading frame structure, and protein-coding potential reasonably well. The "hard" problems for DNA LMs are almost exclusively non-coding.

### 5.2 The Tokenization Problem Is Underappreciated

DNABERT uses k-mer tokenization, creating overlapping tokens and ambiguous representations. DNABERT-2 switched to BPE tokenization. The choice of tokenization dramatically affects what counts as "context" and whether multi-nucleotide motifs fall within a single token or span multiple tokens.

### 5.3 Scale Has Not Solved the Non-Coding Problem

NT-v2 (500M–2.5B parameters) does not dramatically outperform smaller models on regulatory element tasks. Evo2 at 40B parameters shows improvements, but primarily on tasks with strong sequence-level signals. The non-coding regulatory grammar problem appears to be a *data supervision* problem, not a *model capacity* problem. More sequence data does not substitute for functional annotation.

### 5.4 Enformer-Style Supervised Models Often Beat Foundation Models

For most regulatory prediction tasks, Enformer and successors — trained with rich supervised labels — outperform self-supervised DNA LMs when those LMs are fine-tuned. This suggests that for regulatory genomics specifically, the self-supervised pretraining signal is not as informative as supervised training on functional data.

**Implication:** Incorporating scRNA-seq as a supervision signal may be more valuable than scaling the DNA LM further. The researcher's intuition appears well-grounded.

### 5.5 Cross-Species Generalization Is Inconsistent

For non-coding regulatory sequences, where conservation is weaker and regulatory programs are lineage-specific, cross-species training data may actually hurt by providing misleading signal. The regulatory vocabulary is not universal across species in the way that codon usage is.

### 5.6 Evaluation Leakage Concerns

Several DNA LM papers have been criticized for using evaluation sets not cleanly separated from training data by chromosomal holdout. Models with very large context windows are harder to evaluate cleanly.

### 5.7 RNA-Level Supervision Is Not Just About Expression Levels

scRNA-seq provides:
- **Splicing ratios** at single-cell resolution
- **RNA velocity** (nascent vs. spliced RNA via intronic reads) — a direct readout of transcriptional activity
- **Allele-specific expression** in heterozygous individuals — directly links specific sequence variants to expression in specific cells
- **Cellular context for variant effects** — the same variant can have different effects in different cell types

---

## 6. Summary Table of Key Gaps

| Limitation | Severity | Potentially fixable with scRNA-seq? |
|---|---|---|
| Promoter-enhancer pairing | High | Yes (scATAC+scRNA co-variation) |
| Cell-type-specific regulatory activity | High | Yes (direct supervision) |
| Long-range context utilization | High | Partial (training signal; not sufficient alone) |
| Gene-gene regulatory relationships | High | Partial (co-expression) |
| Variant effects in non-coding regions | Medium | Yes (cell-type-specific eQTL signal) |
| Isoform/splicing specificity | Medium | Yes (with long-read scRNA) |
| 3D genome organization | Medium | No (requires Hi-C / structural data) |
| Cross-species non-coding transfer | Medium | No (organism-specific) |
| Repeat/TE-derived regulatory elements | Medium | No (training data issue) |
| Protein domain / gene function | Low (out of scope for DNA LM) | No (requires protein data) |

---

*Key papers to verify: Dalla-Favera et al. (NT-v2), Nguyen et al. (Evo2 Arc Institute), Schiff et al. (Caduceus), Marin et al. (BEND), Avsec et al. (Enformer), Linder et al. (Borzoi), Grešová et al. (Genomic Benchmarks).*
