# GENE_UNDERSTANDING.md

## Do DNA Language Models Learn Gene-Level Concepts and Inter-Gene Relationships?

*A survey for scRNA-seq × DNA LM integration research at Arc Institute*
*Compiled from training knowledge (cutoff August 2025). Specific numerical claims are flagged where memory confidence is moderate.*

---

## 1. Mechanistic Interpretability: Do Attention Heads Correspond to Gene Structures?

### What has been studied

The mechanistic interpretability literature for DNA language models is thin compared to protein LMs (e.g., ESM-2 probing work) and NLP transformers, but a few concrete findings exist.

**Nucleotide Transformer (NTv3 family — Dalla-Favera et al. / InstaDeep + EMBL-EBI, 2023)**
The original Nucleotide Transformer paper evaluated representations via downstream task performance rather than mechanistic dissection. They showed that embeddings from intermediate layers outperform final-layer embeddings for regulatory element prediction tasks (promoter, enhancer, splice site classification). This "middle-layer peak" pattern is common across biological sequence LMs and suggests the model builds up structural knowledge in early-to-mid layers before collapsing it toward a next-token prediction objective in later layers. The paper does not contain attention-head–level dissection.

**Hyena / HyenaDNA (Nguyen et al., 2023, *NeurIPS*)**
HyenaDNA's authors demonstrated that models trained with context lengths up to 1 Mb can encode positional relationships spanning gene-scale distances. Single-nucleus chromatin accessibility prediction improved substantially at longer context, implying the model captures distal regulatory relationships. However, this is a behavioral result, not a mechanistic one — there is no identification of which convolutional filters or layers correspond to specific gene features.

**Caduceus / MambaDNA**
State-space model architectures (Caduceus, 2024) focus on benchmark performance; no attention-map analysis is possible because these architectures lack attention heads.

**Evo2 (Arc Institute / Together AI, 2024–2025)**
Evo2 was trained on ~9 trillion nucleotides across all domains of life. The preprint describes generating novel functional proteins and predicting variant effects, but does not contain mechanistic interpretability analyses mapping specific layers to biological features like TSS position or exon boundaries. Given Evo2's use of a StripedHyena architecture, classical attention-head probing is architecturally inapplicable.

**DNABERT family**
DNABERT (Ji et al., 2021, *Bioinformatics*) applied k-mer tokenization to human reference genome segments. The paper includes attention visualization showing that attention weights cluster around known regulatory motifs (TATA box, GC-rich promoter regions). This is the closest thing to published "attention head = gene feature" analysis in the DNA LM literature — but it is qualitative visualization, not quantitative circuit-level analysis. DNABERT-2 (Zhou et al., 2023) replaced k-mer tokenization with BPE and showed improved performance but similarly lacks rigorous mechanistic analysis.

**The honest summary:** No paper as of August 2025 has performed the level of mechanistic interpretability on a DNA LM that has been done for ESM-2 (which has papers identifying attention heads corresponding to contact maps) or GPT-2 (induction heads, etc.). The field has relied on probing classifiers and benchmark tasks as proxies.

---

## 2. Probing Studies: Gene Identity, Expression Level, Regulatory State

### Gene identity

**Does a model know it is "inside gene X"?**
There are no published probing studies asking a DNA LM to recover gene identity (i.e., Entrez gene ID or gene name) from internal representations. This is a meaningful gap. Gene identity is not a property of local sequence — it requires knowing global position, which models with limited context windows cannot reliably encode. Longer-context models (HyenaDNA at 1 Mb, Evo2 at up to 131 kb in practice) could in principle encode this, but the probing experiment has not been published.

### Expression level

**Can DNA LMs predict RNA expression levels?**

- DNA LMs trained with a pure MLM or next-token-prediction objective do **not** predict steady-state expression levels well in a zero-shot fashion. Expression is determined by chromatin state, trans-acting factors, and cell-type context — not purely by DNA sequence.

- **Enformer (Avsec et al., 2021, *Nature Methods*)** is the clearest counter-example, but Enformer is not a self-supervised LM. It is a supervised sequence-to-function model trained on CAGE, ATAC, ChIP-seq, and RNA-seq tracks from ENCODE/Roadmap. It explicitly learns to predict expression from 200 kb of surrounding sequence. Enformer representations *do* contain expression-level information — that is what they were trained on. Probing Enformer's penultimate layer for expression recovers excellent predictions (R² ≈ 0.7–0.85 for many gene–cell-type pairs, from memory; verify in paper).

- **Borzoi (Linder et al., 2023/2024, *Nature Genetics*)** extends Enformer to longer context (524 kb) and RNA-seq at nucleotide resolution. It can predict splicing patterns in addition to expression levels, and probing its representations recovers both.

- For pure DNA LMs (NTv3, Evo2), published variant effect scores correlate with pathogenic variants (ClinVar, gnomAD) but not straightforwardly with expression levels per se.

### Regulatory state / epigenetic state

- **Nucleotide Transformer (2023):** Downstream probing on 18 tasks including histone modification prediction (H3K4me3, H3K27ac, etc.) and open chromatin (ATAC). Mid-layer representations achieve strong accuracy. This confirms the model has learned features correlated with regulatory state — but this is correlation, not mechanism.

- **GENA-LM (Fishman et al., 2023, *bioRxiv*):** Systematic probing on a broader task set including promoter recognition, splice site prediction, and enhancer classification. Similar "mid-layer peak" pattern.

- **Important caveat:** Regulatory state probing tells us the representations *contain information about* regulatory features, not that the model has an internal representation of "this region is an active enhancer in liver." The model may have memorized sequence motifs that are statistically associated with these labels.

---

## 3. Co-expression and Gene Network Learning

### The core hypothesis and why it is almost certainly correct

Your hypothesis — that pure DNA LMs **cannot** learn co-expression relationships or regulatory networks — is well-founded for the following reasons:

1. **Co-expression is a trans-cellular-state property.** Two genes are co-expressed because in the same cell, under the same conditions, their transcripts appear together. This information is simply not present in a DNA sequence training corpus.

2. **Gene regulatory networks require perturbation data** (or at minimum, multi-condition expression data). A model trained on sequence cannot observe that knocking out TF X reduces expression of gene Y.

3. **Synteny is not co-expression.** Genomic proximity (useful for predicting cis-regulatory relationships) does not predict co-expression. Functionally related genes are often not syntenic in mammals.

### What DNA LMs *do* capture that is adjacent to this

- **Shared regulatory motifs:** If genes A and B both have CTCF binding sites and AP-1 motifs in their promoters, a DNA LM may assign similar embeddings to those promoter regions. This is a very weak proxy for co-regulation — it captures TF co-binding but not condition-specific co-expression.

- **Gene family relationships:** Paralogs share sequence, so DNA LM embeddings of paralogous gene bodies will cluster together. This is sequence homology, not regulatory learning.

- **TAD-scale context:** Long-context models (HyenaDNA, Evo2) can process entire TADs, which are the scale at which co-regulated gene clusters exist. There is no published evidence that this leads to learning co-expression, but it is a theoretical prerequisite. A probing study specifically asking "do Evo2 representations of gene A predict the expression of gene B in the same TAD" would be a meaningful experiment and has not been done.

### Published negative results / null comparisons

No paper explicitly tests a DNA LM on a co-expression prediction task and reports negative results. This is itself a finding: **the evaluation landscape has not asked this question.** Papers evaluate DNA LMs on sequence-label tasks (variant effects, regulatory element classification) — not on relational tasks between genes.

---

## 4. Expression-Supervised Models vs. Pure DNA LMs

### Enformer (Avsec et al., 2021)

- **Training objective:** Predict 5,313 genomic tracks (CAGE, histone marks, ATAC, RNA-seq, etc.) from 200 kb sequence windows, supervised on ENCODE/Roadmap data across ~100 human cell types + mouse.
- **What it learns:** Direct cell-type–specific expression prediction. TF binding sites, enhancer-promoter links, splice site usage.
- **Inter-gene relationships:** Enformer can predict expression of a gene given the sequence context of its neighbors (it sees ~100 kb upstream and downstream). This allows it to implicitly learn that nearby genes sharing regulatory elements are co-regulated — but only within the 200 kb window.
- **Limitation:** Does not model trans-regulation, does not generalize to cell types outside its training distribution without fine-tuning.

### Borzoi (Linder et al., 2023/2024)

- Extends Enformer to RNA-seq at nucleotide resolution. The internal representations contain information about transcript architecture, not just promoter activity.
- **Published probing result:** ISM on Borzoi identifies known splice regulatory elements (ESEs, ESSs) without any explicit splice annotation in training — the model learned them from RNA-seq signal alone. This is a meaningful "emergent gene structure learning" result.

### scFoundation / Geneformer / scGPT (expression-data LMs, for contrast)

These models are trained on scRNA-seq data directly and are the most direct comparison point:

- **Geneformer (Theodoris et al., 2023, *Nature*):** Trained on ~30 million single-cell transcriptomes. Rank-based tokenization of gene expression. Fine-tuned Geneformer correctly predicts dosage-sensitive genes and, via in silico perturbation of TF nodes, produces biologically coherent downstream expression changes. This is the clearest published evidence that expression-supervised models learn gene networks that DNA LMs cannot.

- **scGPT (Cui et al., 2024, *Nature Methods*):** Similarly trained on millions of cells. Learns cell-type–specific gene programs and perturbation responses.

- **The key contrast:** Geneformer and scGPT tokens are *genes*, and their context window is a cell's entire transcriptome. The model literally observes which genes are expressed together across millions of cells. A DNA LM's tokens are nucleotides, and its context window is a genomic interval. These are fundamentally different learning signals.

### Tabular summary (from memory — verify specific numbers in papers)

| Model | Training signal | Learns expression? | Learns co-expression? | Learns gene regulatory networks? |
|---|---|---|---|---|
| NTv3 / DNABERT | DNA sequence (MLM/NTP) | Weakly, via proxy tasks | No | No |
| Evo2 | DNA sequence (NTP, multi-kingdom) | Variant effects, not expression levels | No evidence | No |
| HyenaDNA | DNA sequence (NTP) | No | No | No |
| Enformer | Supervised CAGE/ChIP/ATAC | Yes, cell-type specific | Within cis window only | Cis-regulatory only |
| Borzoi | Supervised RNA-seq / CAGE | Yes + splicing | Within cis window only | Cis-regulatory only |
| Geneformer | scRNA-seq transcriptomes | Yes (trained on it) | Yes | Yes, via in silico perturbation |
| scGPT | scRNA-seq transcriptomes | Yes (trained on it) | Yes | Partially |

---

## 5. Surprising Findings and Gaps in the Evaluation Landscape

### 5.1 No "gene-level" probing benchmark exists for DNA LMs

The Genomics Benchmark and BEND benchmark evaluate regulatory element classification, variant effect scoring, and splice site prediction — all *local sequence* tasks. No benchmark asks: "does the model's representation of gene X differ systematically from gene Y, in a way that captures gene function?" This is a real gap.

### 5.2 The "gene body" as a unit is not privileged in DNA LM training

These models see a stream of nucleotides. Gene boundaries are not marked. A 512-token window might contain half an intron. There is no reason the model would develop a representation of "a gene" as a coherent unit unless gene-scale structure is necessary to predict local sequence.

### 5.3 Evo2's evolutionary conservation signal may be the closest analog to gene learning

Evo2 was trained across prokaryotes, eukaryotes, viruses, and phage. Conserved genes should have more consistent representations across species than non-conserved sequence. This means Evo2's gene representations are organized around evolutionary function — not expression programs. These partially overlap but are not the same concept.

### 5.4 The "genomic LM + expression data" fusion has barely been explored

The natural experiment — taking a pretrained DNA LM and probing it for co-expression information before and after fine-tuning on scRNA-seq — has not, to my knowledge, been published as of mid-2025. Papers like **GET (Ge et al., 2024, preprint)** take steps in this direction by using Enformer-style representations as inputs to models trained on single-cell data, but the mechanistic question of whether the DNA LM component "learns" co-expression through the fine-tuning signal has not been dissected.

### 5.5 Regulatory grammar vs. expression state

A DNA LM may learn *cis-regulatory grammar* (motif combinations, spacing, orientation rules) without learning *expression states*. These are distinct. A model can know that "TATA box + initiator + downstream promoter element = strong promoter" without knowing whether that promoter is active in any given cell. Enformer collapses this distinction but pure DNA LMs cannot.

### 5.6 Long-range regulatory relationships are theoretically accessible but practically understudied

Super-enhancers driving oncogenes can be hundreds of kilobases from the TSS. Evo2 and HyenaDNA can in principle see such distances. Whether long-range chromatin looping relationships are encoded in the sequence itself (via CTCF motif pairing, cohesin loading sites, etc.) is an open question. Some evidence from CTCF motif grammar studies (Orca, C.Origami) suggests yes, partially — but this has not been connected back to DNA LM representations.

---

## Summary for Research Direction

The evidence as of mid-2025 supports the researcher's hypothesis strongly:

1. **Pure DNA LMs (NTv3, Evo2, HyenaDNA) do not learn co-expression or gene regulatory networks.** The training signal simply does not contain this information.

2. **DNA LMs do learn cis-regulatory features** (promoter motifs, splice sites, conservation) that are partially predictive of expression, but only in the sense that they capture sequence grammar — not expression state.

3. **Expression-supervised models (Enformer, Borzoi) learn cell-type–specific cis-regulatory relationships** but are still fundamentally limited to the cis window and cannot capture trans-regulation.

4. **Expression-data LMs (Geneformer, scGPT) are the only class that demonstrably learns co-expression and gene networks,** and they do so because they see gene co-occurrence across millions of cells.

5. **The most interesting unexplored territory** is whether a DNA LM fine-tuned or jointly trained with scRNA-seq data acquires expression-level representations — and whether such a model's sequence representations become better calibrated to gene regulatory structure.

---

## Key References

- Dalla-Favera et al. (2023). "The Nucleotide Transformer: Building and Evaluating Robust Foundation Models for Human Genomics." *Nature Methods*.
- Ji et al. (2021). "DNABERT: pre-trained Bidirectional Encoder Representations from Transformers model for DNA-language in genome." *Bioinformatics*.
- Zhou et al. (2023). "DNABERT-2: Efficient Foundation Model and Benchmark For Multi-Species Genome." *arXiv:2306.15006*.
- Nguyen et al. (2023). "HyenaDNA: Long-Range Genomic Sequence Model at Single Nucleotide Resolution." *NeurIPS 2023*.
- Avsec et al. (2021). "Effective gene expression prediction from sequence by integrating long-range interactions." *Nature Methods*.
- Linder et al. (2023/2024). "Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation." *Nature Genetics* (verify exact venue).
- Theodoris et al. (2023). "Transfer learning enables predictions in network biology." *Nature*.
- Cui et al. (2024). "scGPT: toward building a foundation model for single-cell multi-omics using generative AI." *Nature Methods*.
- Greššová et al. (2023). "Genomic benchmarks: a collection of datasets for genomic sequence classification." *BMC Genomic Data*.
- Marin et al. (2023). "BEND: Benchmarking DNA Language Models on biologically meaningful tasks." *arXiv:2311.12570*.
- Nguyen et al. (2024). "Sequence modeling and design from molecular to genome scale with Evo." *Science* (Evo1; Evo2 preprint 2025).

*Note: Specific R² values, benchmark numbers, and some publication venues are from training memory and should be verified against the primary papers before citing.*
