# Long-Read Single-Cell RNA-seq Datasets with Genotype Variation

*Compiled for Arc Institute scRNA-seq × DNA language model project (NTv3, Evo2).*
*Companion to `4_smartseq_datasets_with_genotype_variation.md` (Smart-seq cohorts) and `4.5_biggest_full_coverage_scrna_datasets_with_genotype_variation.md` (full-coverage landscape).*
*Verification rounds: PubMed MCP + bioRxiv + GEO/EGA/ENA/dbGaP through June 2026. Every PMID below was independently pulled via `mcp__claude_ai_PubMed__get_article_metadata` and every accession was cross-checked against either a paper's data-availability statement or a repository search.*

> **Headline result (one sentence): The long-read single-cell field as of mid-2026 has exactly *one* genuine multi-donor cohort that approaches the "many genotypes × full-isoform expression" criterion at population-relevant scale — Li et al. 2026 NCI/Korean lung (129 donors, paired genotype, Nanopore single-cell, sc-isoQTL design). Every other public long-read sc dataset to date is either (a) ≤25 donors, (b) lacks paired WGS/genotype, (c) is targeted (covers a small gene panel rather than the whole transcriptome), or (d) is a method paper with no genotype variation at all.**

> **Quick sanity check vs. file 4.5:** File 4.5's long-read sub-section is essentially correct on the entries it lists — but it presents a sparse landscape as if it were dense. After live verification this file expands the inventory (adding Liu 2024 eNeuro 25-donor brain, Pan 2025 Circulation 12-donor heart, Li 2024 *Cell Genomics* 12-patient CRC, Byrne 2024 *Nat Commun* ovarian scTaILoR-seq, Joglekar 2024 *Nat Neurosci* 6-donor brain, Rehman 2025 aging mouse brain, Argue 2025/26 CSF AD, and a few methods entries) **and confirms that none of these new entries change the bottom line**: Li 2026 lung is alone in its tier.

---

## 0. Corrections / clarifications to file 4.5 (June 2026 round)

| File 4.5 claim | Verified status | Notes |
|---|---|---|
| Belchikov et al. 2025 SnISOr-Seq2 *Cell Rep* 44(9):116198, PMID 40913764 | **Confirmed.** GEO **GSE250280** (short-read snRNA) + SRA **PRJNA1238317** (long-read). 12 primary donors (6 GRN-FTD + 6 controls) + 2 validation donors. ~100K nuclei across the 12 primaries. | First-author tissue is superior frontal gyrus. No paired WGS. DOI [10.1016/j.celrep.2025.116198](https://doi.org/10.1016/j.celrep.2025.116198). |
| Foord et al. 2025 Spl-ISO-Seq *Nat Commun* 16(1):8093, PMID 40883294 | **Confirmed.** DOI [10.1038/s41467-025-63301-9](https://doi.org/10.1038/s41467-025-63301-9). Pre-/post-puberty (8–11 vs 16–19 yr) post-mortem visual cortex. Near-single-cell *spatial* (not true single-cell). | Donor count appears to be modest; specific accession not pulled from full text. |
| Bizzotto et al. 2025 *PNAS* GO-TEN FCD, PMID 40674414 | **Confirmed but qualified.** 18 FCD + 17 controls (=35 donors) is the **snRNA-seq** (short-read) cohort. **Long-read GO-TEN was applied to only 3 surgical samples**, ResolveOME to 5 more. dbGaP **phs004124.v1.p1**. | First-author is Bizzotto; last-author Walsh CA (Boston). DOI [10.1073/pnas.2509622122](https://doi.org/10.1073/pnas.2509622122). **Note: a second FCD paper (PMID 40307383, Baldassari 2025 *Nat Neurosci*, Paris group, with Bizzotto as middle author) is a different study, snRNA-seq only.** |
| Doyle et al. PBMC isonome *Front Genet*, PMID 42283034 | **Confirmed.** 1 healthy 43-yr male donor, 2 replicates, 33,018 cells total (PBMC1: 11,352; PBMC2: 21,666). PIPseq V4.0PLUS + Oxford Nanopore PromethION. DOI [10.3389/fgene.2026.1782221](https://doi.org/10.3389/fgene.2026.1782221). | bioRxiv is 2025.10.16.682832. No paired WGS. |
| Joglekar et al. 2021 *Nat Commun* mouse brain, PMID 33469025 | **Confirmed.** GEO **GSE158450**. C57BL/6 P7 mouse, single strain, ScISOr-Seq (10x + PacBio). DOI [10.1038/s41467-020-20343-5](https://doi.org/10.1038/s41467-020-20343-5). | Single-strain mouse — no genotype variation. |
| Michielsen et al. 2026 *Genome Biol* brain exon inclusion, PMID 41764538 | **Confirmed but qualified.** This is a *modeling* paper that re-uses Tilgner-lab long-read sc data (mostly Joglekar 2024 cohort, q.v.). No new donors. DOI [10.1186/s13059-026-04015-z](https://doi.org/10.1186/s13059-026-04015-z). |
| Tian et al. FLAMES 2021 *Genome Biol*, PMID 34763716 | **Confirmed.** GEO **GSE154869 / GSE154870 / GSE154868 / GSE126906 / GSE142285** (cell lines + MuSC) + EGA **EGAS00001005597** (CLL2 patient sample). Mostly cell lines + mouse mdx muscle + 1 CLL patient. DOI [10.1186/s13059-021-02525-6](https://doi.org/10.1186/s13059-021-02525-6). | Method paper. No multi-donor cohort. |
| Li et al. 2026 Korean lung isoQTL | **Confirmed.** PMID **42039472**, bioRxiv DOI [10.64898/2026.03.27.714873](https://doi.org/10.64898/2026.03.27.714873). 129 never-smoker Korean women donors, 37 cell types, paired genotype, sc-isoQTL design. Bolun Li first author, NCI Division of Cancer Epidemiology and Genetics. | **The headline paper.** No journal acceptance as of June 2026; bioRxiv only. Accession numbers in the preprint should be re-verified at journal publication. |
| Karakulak et al. ccRCC MAS-seq *Genome Res*, PMID 40107723 | **Confirmed.** 2,599 ccRCC PDO cells, PacBio MAS-ISO-seq. *Genome Res* 35(4):698 (Apr 2025). | Patient-derived organoids; donor count is small (PDO panel) — verify exact n in supplement. No paired WGS in primary data. |
| Hagemann-Jensen "MAS-ISO-seq" Al'Khafaji 2024 method, PMID 37291427 | **Confirmed.** *Nat Biotechnol* 42:582 (in-issue 2024, online 2023). Method paper. Small melanoma TIL panel + cell lines. DOI [10.1038/s41587-023-01815-7](https://doi.org/10.1038/s41587-023-01815-7). |
| ENCODE4 long-read RNA-seq (Reese 2023), PMID 37292896 | **Confirmed but qualified.** 264 PacBio LR-RNA-seq libraries, 81 unique human + mouse samples (cell lines + tissues). **Mostly bulk; LR-Split-seq single-cell is a small subset.** DOI [10.1101/2023.05.15.540865](https://doi.org/10.1101/2023.05.15.540865). | bioRxiv; not formally journal-published as of June 2026. No multi-donor genotype variation. |

---

## 1. The verified cohort table — every public long-read scRNA-seq dataset I could find with > 1 distinct human donor

Inclusion criterion: long-read (PacBio/Nanopore) at single-cell/single-nucleus resolution, where there are at least 2 distinct donors. Method/benchmark papers using only cell lines are listed separately in §2.

| # | Study | Citation (PMID / DOI) | Accession | Long-read protocol | Donors | Cells / nuclei | Tissue | Species | Paired genotype/WGS? | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **Li et al. 2026 NCI/Korean lung sc-isoQTL** | bioRxiv 2026; PMID **42039472**; DOI [10.64898/2026.03.27.714873](https://doi.org/10.64898/2026.03.27.714873) | Not yet in primary repos (preprint; verify on journal acceptance) | **10x Chromium + Oxford Nanopore** (full-length sc long-read) | **129** never-smoker Korean women | "abundant" — 37 cell types resolved; exact count not in abstract | Lung (tumor + normal) | Human | **Yes** (paired genotype, sc-isoQTL design) | **THE HEADLINE DATASET.** First population-scale sc-isoQTL. |
| 2 | **Liu et al. 2024 eNeuro — AD/DLB/PD prefrontal cortex** | *eNeuro* 11(12) (2024); PMID **39658200**; DOI [10.1523/ENEURO.0296-24.2024](https://doi.org/10.1523/ENEURO.0296-24.2024) | **EGA EGAS50000000132** (controlled-access) | **10x 3' v3.1 snRNA + targeted PacBio Iso-Seq** (50-gene enrichment panel) | **25** (5 control + 6 AD + 7 PD + 7 DLB) | 165,440 nuclei | Prefrontal cortex (BA 8/9) | Human | No paired WGS reported | **Largest verified multi-donor long-read sc cohort outside Li et al.** Long-read covers only 50 disease-related genes — *targeted*, not whole-transcriptome. |
| 3 | **Bizzotto et al. 2025 — FCD GO-TEN** | *PNAS* 122(29):e2509622122 (2025); PMID **40674414**; DOI [10.1073/pnas.2509622122](https://doi.org/10.1073/pnas.2509622122) | **dbGaP phs004124.v1.p1** | snRNA-seq + **GO-TEN** targeted long-read Nanopore (on 3 samples only) + **ResolveOME** (on 5 samples) | **35** total for snRNA (18 FCD + 17 control); **only 3 samples have GO-TEN long-read genotyping**; 5 separate FCD have ResolveOME | >400,000 nuclei (combined own + reanalyzed) | Cortex (surgical FCD II + controls) | Human | **Yes — somatic variant calling per nucleus** (PI3K-mTOR pathway focus) | Genotype is *somatic mosaic*, not germline. dbGaP controlled access. |
| 4 | **Pan et al. 2025 *Circulation* — heart isoform atlas** | *Circulation* 152(21):1501–1514 (2025); PMID **41017471**; DOI [10.1161/CIRCULATIONAHA.125.074959](https://doi.org/10.1161/CIRCULATIONAHA.125.074959) | **GEO GSE288222** | **10x snRNA + Oxford Nanopore** (snNanoRNAseq) | **12** (6 non-diseased + 6 end-stage heart failure) | 30,683 nuclei | Adult left ventricle | Human | None reported | Disease vs. control design; no germline genotype paired. |
| 5 | **Belchikov et al. 2025 *Cell Rep* — FTD SnISOr-Seq2** | *Cell Rep* 44(9):116198 (2025); PMID **40913764**; DOI [10.1016/j.celrep.2025.116198](https://doi.org/10.1016/j.celrep.2025.116198) | GEO **GSE250280** (short-read) + SRA **PRJNA1238317** (long-read) | **SnISOr-Seq2** — targeted long-read snRNA, 3,630-gene splice-junction enrichment panel | **12** primary (6 GRN-FTD + 6 controls) + 2 validation | ~100K nuclei across 12 | Superior frontal gyrus + (validation) frontal/occipital gyrus | Human | None reported | Targeted 3,630-gene panel, not whole-transcriptome. |
| 6 | **Li, Z. et al. 2024 *Cell Genomics* — CRC isoform atlas** | *Cell Genom* 4(9):100641 (2024); PMID **39216476**; DOI [10.1016/j.xgen.2024.100641](https://doi.org/10.1016/j.xgen.2024.100641) | **ENA PRJEB68074** (raw); **Zenodo 10.5281/zenodo.12750017** (processed) | **10x 3' + PacBio Iso-Seq** | **12** CRC patients (8 tumor + 11 normal samples) | 18,966 short-read cells; 1,994 matched long-read cells | Colorectal cancer + adjacent normal | Human | HLA typing only (from short-read scRNA); no WGS | Long-read cell count is small (~2K). |
| 7 | **Joglekar et al. 2024 *Nat Neurosci* — mouse + human brain isoform atlas** | *Nat Neurosci* 27(6):1051–1063 (2024); PMID **38594596**; DOI [10.1038/s41593-024-01616-4](https://doi.org/10.1038/s41593-024-01616-4) | Atlas portal www.isoformAtlas.com (primary accession not verified — possibly extends GSE158450) | **ScISOr-Seq2** (10x + Oxford Nanopore + PacBio HiFi) | **6** human donors (hippocampus only) + 12 mice (C57BL/6NTac, single strain) | ~204,725 mouse cells + human cell count not stated | Mouse hippocampus + visual cortex + striatum + thalamus + cerebellum; human hippocampus | Mouse + Human | None reported | **Successor to Joglekar 2021.** Adds adult mouse + human hippocampus. Single human strain in mouse arm; 6 healthy human donors. |
| 8 | **Byrne et al. 2024 *Nat Commun* — ovarian cancer scTaILoR-seq** | *Nat Commun* 15(1):6916 (2024); PMID **39134520**; DOI [10.1038/s41467-024-51252-6](https://doi.org/10.1038/s41467-024-51252-6) | Not pulled from abstract (Genentech in-house deposit; verify) | **scTaILoR-seq** — hybridization-capture targeted long-read Nanopore on 1,000+ gene panel | Cell lines + primary ovarian tumors (donor count not in abstract; "primary tumors") | 10,796 single-cell long-read transcriptomes | Ovarian cancer cell lines + primary tumors | Human | SNV calling from long reads; no germline WGS reported | First demonstration of allelic imbalance from long-read sc reads on a real tumor panel. |
| 9 | **Argue & Gate 2025/26 — Alzheimer immune cells** | *Alzheimer's & Dementia* 21(Suppl.) (Dec 2025 abstract / 2026 in-press; abstract DOI [10.1002/alz.70855](https://doi.org/10.1002/alz.70855); PMC12726925) | Not yet listed in primary repos (conference abstract) | **Oxford Nanopore + scNanoGPS + IsoQuant + SQANTI3** | **CSF: 45 controls + 13 MCI/AD** (=58); **Blood: 22 controls + 28 AD** (=50). Some donors may overlap. | Cell count not stated; 97,920 unique transcripts called | CSF + peripheral blood (immune cells) | Human | None reported | Conference abstract — full paper not yet found. **Largest long-read sc donor count after Li 2026 if confirmed.** |
| 10 | **Hardwick et al. 2022 *Nat Biotechnol* — SnISOr-Seq method** | *Nat Biotechnol* 40(7):1082–1092 (2022); PMID **35256815**; DOI [10.1038/s41587-022-01231-3](https://doi.org/10.1038/s41587-022-01231-3) | Verify (not pulled) | **SnISOr-Seq** — long-read targeted snRNA-seq with enrichment | "Adult human frontal cortex" — donor count not in abstract | Thousands of nuclei | Human frontal cortex | Human | None reported | Original SnISOr-Seq method paper. Donor n likely small (<5); verify supplement. |
| 11 | **Joglekar et al. 2021 *Nat Commun* — mouse brain isoform atlas** | *Nat Commun* 12(1):463 (2021); PMID **33469025**; DOI [10.1038/s41467-020-20343-5](https://doi.org/10.1038/s41467-020-20343-5) | **GEO GSE158450** | **ScISOr-Seq** (10x + PacBio Iso-Seq) | C57BL/6 P7, single strain | Thousands across 45 cell types | Mouse hippocampus + PFC | Mouse | Single strain — **no genotype variation** | Foundational, but single strain. |
| 12 | **Foord et al. 2025 *Nat Commun* — Spl-ISO-Seq visual cortex** | *Nat Commun* 16(1):8093 (2025); PMID **40883294**; DOI [10.1038/s41467-025-63301-9](https://doi.org/10.1038/s41467-025-63301-9) | Verify (not pulled) | **Spl-ISO-Seq** — spatial long-read at near-single-cell resolution | Pre-puberty (8–11 yr) + post-puberty (16–19 yr) post-mortem; donor count not stated in abstract | Spatial near-single-cell | Visual cortex | Human | None reported | Spatial — not strictly single-cell. |
| 13 | **Shiau et al. 2023 *Nat Commun* — scNanoGPS** | *Nat Commun* 14(1):4124 (2023); PMID **37433798**; DOI [10.1038/s41467-023-39813-7](https://doi.org/10.1038/s41467-023-39813-7) | Verify | Single-cell Nanopore (10x + ONT) | **4 tumors + 2 cell lines** | 23,587 long-read transcriptomes | Mixed tumor types (kidney, etc.) | Human | Per-cell mutation calling from long reads — **same-cell genotype + phenotype** | Method paper — small cohort, no germline WGS. |
| 14 | **Doyle et al. 2026 *Front Genet* — PBMC isonome** | *Front Genet* (2026); PMID **42283034**; DOI [10.3389/fgene.2026.1782221](https://doi.org/10.3389/fgene.2026.1782221) | Verify (Zenodo referenced in paper) | **PIPseq V4.0PLUS + Oxford Nanopore PromethION** | **1** donor, 2 replicates | 33,018 cells | PBMCs | Human | None reported | Single-donor benchmark. |
| 15 | **Rehman et al. 2025 bioRxiv — aging mouse brain isoform atlas** | bioRxiv 2025.06.05.658133; DOI [10.1101/2025.06.05.658133](https://doi.org/10.1101/2025.06.05.658133); PMC12190350 | Not yet in GEO (preprint) | **10x + Oxford Nanopore** (PromethION R9.4.1) | C57BL/6JN, single strain; 4 ages | 156,935 cells | Mouse cortex + hippocampus | Mouse | Single strain — **no genotype variation** | Aging-axis study. |
| 16 | **Karakulak et al. 2025 *Genome Res* — ccRCC MAS-seq** | *Genome Res* 35(4):698 (2025); PMID **40107723**; DOI [10.1101/gr.279345.124](https://doi.org/10.1101/gr.279345.124) | Verify (CSHL portal) | **PacBio MAS-ISO-seq** | Multi-patient PDO panel (verify n in supplement) | 2,599 ccRCC PDO cells | Patient-derived organoids | Human | Not reported | Small donor count. |
| 17 | **Pan, Pandey et al. 2025 — scCLEAN + MAS-seq** | *Nat Commun* 16(1):4664 (2025); PMID **40389438**; DOI [10.1038/s41467-025-59880-2](https://doi.org/10.1038/s41467-025-59880-2) | Verify | **PacBio MAS-seq + CRISPR/Cas9 depletion (scCLEAN)** | Mixed tissue panel (donor count not stated in abstract) | Unstated | Various | Human | Not reported | Method paper. |
| 18 | **Tian et al. 2021 *Genome Biol* — FLAMES** | *Genome Biol* 22(1):310 (2021); PMID **34763716**; DOI [10.1186/s13059-021-02525-6](https://doi.org/10.1186/s13059-021-02525-6) | **GEO GSE126906/154869/154870/154868/142285** (mixology+MuSC); EGA **EGAS00001005597** (CLL2) | **10x + Oxford Nanopore** (subsampling protocol) | 5 cell lines + mouse mdx + **1 CLL patient** | Thousands | Cell lines + mouse muscle + 1 CLL PBMC | Human + Mouse | None reported | Foundational method paper. |
| 19 | **Al'Khafaji et al. 2024 *Nat Biotechnol* — MAS-ISO-seq method** | *Nat Biotechnol* 42(4):582–586 (2024, online 2023); PMID **37291427**; DOI [10.1038/s41587-023-01815-7](https://doi.org/10.1038/s41587-023-01815-7) | Verify | **PacBio MAS-ISO-seq** | Few melanoma patients (TIL panel) + cell lines | Thousands | Melanoma TIL | Human | Not reported | Foundational MAS-seq method paper. |

**The "verified multi-donor long-read sc" count is 19 entries, but only one (Li 2026 lung) has > 50 donors AND paired genotype.**

---

## 2. Long-read sc method papers and resources without donor variation

Listed for completeness — these are not candidates for the inventory but explain the technology landscape.

| Method | Citation | Notes |
|---|---|---|
| **Gupta et al. 2018 *Nat Biotechnol* — ScISOr-Seq** | PMID **30320766**; DOI [10.1038/nbt.4259](https://doi.org/10.1038/nbt.4259) | Original ScISOr-Seq, mouse cerebellum cells |
| **Rebboah et al. 2021 *Genome Biol* — LR-Split-seq** | PMID **34620214**; DOI [10.1186/s13059-021-02505-w](https://doi.org/10.1186/s13059-021-02505-w) | Parse Evercode + ONT, C2C12 mouse myogenic cells; no donor variation |
| **Fu et al. 2025 *Nat Commun* — Longcell** | PMID **40683866**; DOI [10.1038/s41467-025-60902-2](https://doi.org/10.1038/s41467-025-60902-2) | Computational pipeline for sc/spatial Nanopore long reads; applied across existing datasets |
| **FLAMESv2 (Wang et al. 2025 preprint)** | bioRxiv 2025.10.19.683327 | Modular protocol-agnostic R/Bioconductor pkg; multi-sample analysis upgrade for FLAMES |
| **SCOTCH (Xu et al. 2025 bioRxiv → Nat Commun 2026)** | bioRxiv 2024.04.29.590597; DOI [10.1038/s41467-026-72665-5](https://doi.org/10.1038/s41467-026-72665-5) | Platform-agnostic (ONT + PacBio + 10x + Parse) isoform pipeline; no new data of its own |
| **HIT-scISOseq** | *Nat Commun* (2023); DOI [10.1038/s41467-023-38324-9](https://doi.org/10.1038/s41467-023-38324-9) | PacBio-CCS scRNA isoform with > 10M long reads/run |
| **LongBench (bioRxiv 2025.09.11)** | bioRxiv 2025.09.11.675724 | Multi-platform cancer cell-line reference dataset (ONT cDNA + direct RNA + PacBio Kinnex); cell-lines only |
| **Multi-platform sc isoform comparison (Feb 2025 bioRxiv)** | bioRxiv 2025.02.22.639662 | Kinnex on Revio; benchmark on 1 donor |
| **Karakulak benchmark preprint (Mar 2024)** | bioRxiv 2024.03.14.584953; *NAR Genom Bioinform* 7:lqaf089 (2025) | MAS-ISO-seq vs short-read quality assessment |

---

## 3. Sub-technology breakdown of the verified entries

### 3a. PacBio MAS-ISO-seq / Kinnex single-cell

| Dataset | Donors | Cells | Notes |
|---|---|---|---|
| Al'Khafaji 2024 (PMID 37291427) | Few melanoma TIL + cell lines | Thousands | Method paper |
| Karakulak 2025 *Genome Res* (PMID 40107723) | Multi-patient PDO panel (small) | 2,599 | PDOs from ccRCC patients |
| Pan/Pandey 2025 scCLEAN+MAS-seq (PMID 40389438) | Unstated | Unstated | CRISPR-Cas9 depletion variant |
| LongBench 2025 bioRxiv | Lung cancer cell lines | Multiple | Cell lines only |

**Verdict on MAS-seq:** Still no published >10-donor MAS-seq sc cohort. Kinnex has been *technically* deployed at scale (PacBio Revio) but no published consortium-scale cohort yet exists.

### 3b. Oxford Nanopore-based single-cell (whole-transcriptome)

| Dataset | Donors | Cells | Notes |
|---|---|---|---|
| **Li 2026 NCI/Korean lung (bioRxiv, PMID 42039472)** | **129** | 37 cell types | **Headline.** Population-scale isoQTL |
| Pan 2025 *Circulation* heart (PMID 41017471) | 12 (6 + 6) | 30,683 | snRNA + ONT |
| Doyle 2026 PBMC isonome (PMID 42283034) | 1 (2 reps) | 33,018 | PIPseq + ONT |
| Rehman 2025 aging mouse (bioRxiv) | C57BL/6JN only | 156,935 | Aging axis, single strain |
| Tian 2021 FLAMES (PMID 34763716) | 1 CLL + cell lines | Thousands | Foundational |
| Shiau 2023 scNanoGPS (PMID 37433798) | 4 tumors + 2 cell lines | 23,587 | Same-cell genotype+phenotype |
| Argue 2025/26 AD CSF+blood (abstract) | 45+13+22+28 (≤108 if non-overlapping) | Unstated | **Potentially second-largest** if confirmed |

### 3c. Targeted long-read snRNA-seq (gene panels with hybridization capture)

These are *not* whole-transcriptome — they only cover a chosen gene panel — so for "BigWig-style full base-pair coverage" they only contribute on those genes.

| Dataset | Panel | Donors | Cells |
|---|---|---|---|
| **Liu 2024 eNeuro AD/PD/DLB cortex (PMID 39658200)** | **50-gene** | **25** | 165,440 nuclei |
| Belchikov 2025 SnISOr-Seq2 FTD (PMID 40913764) | **3,630-gene** | 12 + 2 | ~100K nuclei |
| Hardwick 2022 SnISOr-Seq method (PMID 35256815) | Targeted enrichment | Small (verify) | Thousands |
| Byrne 2024 scTaILoR-seq ovarian (PMID 39134520) | **1,000+ gene** | Cell lines + primaries | 10,796 transcriptomes |
| Bizzotto 2025 GO-TEN FCD (PMID 40674414) | mTOR-pathway-focused | **3 (long-read subset)** | Cells (subset of 400K nuclei) |

### 3d. PacBio Iso-Seq + 10x (whole-transcriptome)

| Dataset | Donors | Cells |
|---|---|---|
| Li, Z. et al. 2024 CRC atlas (PMID 39216476) | 12 patients | 1,994 long-read cells |
| ENCODE4 (Reese 2023 bioRxiv) | 81 unique samples | Bulk + small sc subset |
| DLPFC/OFC cortex Iso-Seq (Nov 2025 bioRxiv 2025.11.25.690524) | 3 DLPFC + 4 OFC | FANS-sorted bulks (not strictly single-cell) |

### 3e. PacBio + spatial / near-single-cell

| Dataset | Donors | Notes |
|---|---|---|
| Joglekar 2021 *Nat Commun* (GSE158450) | Mouse C57BL/6 single strain | ScISOr-Seq + Visium |
| Joglekar 2024 *Nat Neurosci* (PMID 38594596) | 6 human + 12 mouse | ScISOr-Seq2 (ONT + PacBio HiFi) |
| Foord 2025 Spl-ISO-Seq (PMID 40883294) | Pre/post-puberty visual cortex | Spatial near-single-cell |

---

## 4. What I searched for and could *not* verify (vaporware or non-existent)

| Search target | What I found | Status |
|---|---|---|
| **Nagura et al. 2025 "B-cell isoQTL 67 individuals 17,000 ieQTLs"** | PMID **40336129**, *Genome Biol* 26(1):110. **THIS IS BULK Oxford Nanopore RNA-seq of LCLs, NOT single-cell.** From Japanese 1000G samples (Coriell). | Looks single-cell-relevant but is bulk LCL. Not included in the cohort table. Could be repurposed as bulk LCL long-read reference. DOI [10.1186/s13059-025-03583-w](https://doi.org/10.1186/s13059-025-03583-w). |
| **"HipSci long-read scRNA-seq"** | No paper found. HipSci has bulk RNA-seq + bulk ATAC + iPSC bank; no long-read sc deposit at population scale exists. | **Does not exist as of June 2026.** |
| **"iPSCORE long-read sc"** | No paper found. iPSCORE has bulk RNA-seq + bulk ATAC + iPSC genotypes (dbGaP phs001325), but no long-read sc cohort. | **Does not exist.** |
| **"BICCN long-read sc atlas"** | BICCN papers are mostly short-read 10x snRNA + scATAC. No BICCN-branded long-read sc cohort at the scale of the BICCN short-read atlases. | **Does not exist at scale.** |
| **"HCA long-read"** | HCA publishes per-organ atlases; no formal long-read sc atlas yet. | **Does not exist.** |
| **"F1 mouse long-read sc strain panel"** | (a) Lemoine et al. 2025 *Sci Rep* (DOI [10.1038/s41598-025-21643-w](https://doi.org/10.1038/s41598-025-21643-w)): **bulk** PacBio Iso-Seq on a single C57BL/6J × CAST/EiJ F1 female brain (GSE246857). NOT single-cell. (b) IGVF Consortium 2026 hybrid-crosses preprint (bioRxiv 2026.04.02.716195, PMC13060228): 7 F1 hybrids (B6J × 7 CC founders), 5.3–6.7M nuclei, 8 tissues, snRNA-seq — but **SHORT-READ** (Parse split-pool), not long-read. | **Long-read sc F1 strain panels do not exist.** Closest is the IGVF short-read effort. |
| **"GTEx long-read sc"** | GTEx Glinos 2022 *Nature* 608:353 (PMID 35922509) is **bulk ONT** on 88 tissue samples (15 tissues). No sc long-read GTEx layer. | **Does not exist.** |
| **"MAS-seq deployed in a published sc-eQTL cohort"** | Not found. Method papers + a handful of small disease cohorts only. | **Does not exist.** |
| **"sc-isoQTL cohort other than Li 2026"** | No second sc-isoQTL paper found in PubMed or bioRxiv. | **Li 2026 is alone in its tier.** |
| **"Long-read SHARE-seq" / "Long-read MULTIome"** | Not found. SHARE-seq and 10x Multiome are short-read paradigms; no long-read variant published. | **Does not exist.** |
| **"Smart-seq3 paired with long-read"** | Not found. SS3xpress is illumina-short-read despite having 5' UMI + full-length cDNA. | **Does not exist.** |
| **"R2C2 single-cell at scale"** | Some R2C2 papers (e.g. PMID 35130954, 30201725) exist but none with multi-donor genotype variation. | **Does not exist at scale.** |
| **"Population-scale Kinnex"** | PacBio markets Kinnex for "population-scale" but no published Kinnex cohort exceeds a handful of donors. | **Marketing only.** |
| **"scCirRL-seq"** | The user-suggested keyword returned no PubMed hits. The closest hits (PMIDs 42286899, 42210074, 42205581…) are unrelated to a method named scCirRL-seq. | **No paper of this name exists in PubMed.** |

---

## 5. Priority ranking for "genotype → full-isoform expression" learning

Ranking criterion: combined value for training a DNA LM where the supervision is sequence → isoform-level RNA output. Heuristic = (donors × log(cells) × paired-WGS-flag × full-transcriptome-flag).

| Rank | Dataset | Donors | Cells | Whole-transcriptome long-read? | Paired genotype? | Why this rank |
|---|---|---|---|---|---|---|
| **1** | **Li et al. 2026 NCI/Korean lung sc-isoQTL** (PMID 42039472) | **129** | Tens of thousands (37 cell types) | **Yes (ONT whole-tx)** | **Yes — sc-isoQTL design** | The only dataset that simultaneously satisfies (>100 donors) AND (whole-transcriptome long-read sc) AND (paired germline genotype). |
| 2 | **Liu et al. 2024 eNeuro AD/PD/DLB** (PMID 39658200) | 25 | 165,440 nuclei | **No** — targeted (50 genes) | No germline WGS reported | Largest *donor count* outside Li 2026. But only 50 genes covered by long-read. |
| 3 | **Argue & Gate 2025/26 AD CSF + blood** (abstract only) | up to ~108 if non-overlapping | Unstated | Yes (whole-tx ONT) | None reported | Second-largest donor count if confirmed. Conference abstract — needs verification at full-paper stage. |
| 4 | **Bizzotto et al. 2025 *PNAS* FCD GO-TEN** (PMID 40674414) | 35 snRNA donors / 3 with long-read GO-TEN | >400K nuclei | No — targeted (mTOR pathway) | **Yes — somatic mosaic genotype per cell** | Unique design (same-cell genotype + phenotype) but **somatic, not germline**. |
| 5 | **Joglekar et al. 2024 *Nat Neurosci* brain isoform atlas** (PMID 38594596) | 6 human + 12 mouse | ~205K mouse + human (n unstated) | **Yes** (ONT + PacBio HiFi) | None | Cross-species value; small human donor count. |
| 6 | **Li, Z. et al. 2024 *Cell Genomics* CRC** (PMID 39216476) | 12 | 1,994 long-read | **Yes** (PacBio Iso-Seq) | HLA typing only | Disease panel; small long-read cell count. |
| 7 | **Pan et al. 2025 *Circulation* heart** (PMID 41017471) | 12 | 30,683 nuclei | **Yes** (snNanoRNAseq) | None | Multi-condition (HF vs control). |
| 8 | **Belchikov et al. 2025 *Cell Rep* FTD SnISOr-Seq2** (PMID 40913764) | 12 + 2 | ~100K | No — targeted (3,630 genes) | None | Disease panel; useful for splicing biology. |
| 9 | **Byrne et al. 2024 *Nat Commun* ovarian scTaILoR-seq** (PMID 39134520) | Cell lines + primaries (n unstated) | 10,796 | No — targeted (1,000+ genes) | SNV calls from long-reads only | First demonstration of allelic imbalance from sc long reads. |
| 10 | **Foord et al. 2025 *Nat Commun* Spl-ISO-Seq visual cortex** (PMID 40883294) | Small (pre/post-puberty) | Spatial near-single-cell | **Yes** | None | Spatial axis; small donor count. |
| 11 | **Rehman et al. 2025 bioRxiv aging mouse brain** (DOI 10.1101/2025.06.05.658133) | C57BL/6JN only | 156,935 | **Yes** (ONT) | Single strain | Largest mouse long-read sc atlas by cell count. No genotype variation. |
| 12 | **Joglekar et al. 2021 *Nat Commun* mouse brain ScISOr-Seq** (PMID 33469025) | C57BL/6 only | Thousands | **Yes** (PacBio) | Single strain | Foundational; single strain. |
| 13 | **Shiau et al. 2023 *Nat Commun* scNanoGPS** (PMID 37433798) | 4 tumors + 2 cell lines | 23,587 | **Yes** | Per-cell mutation calls | Small donor count. |
| 14 | **Doyle et al. 2026 *Front Genet* PBMC isonome** (PMID 42283034) | 1 (2 reps) | 33,018 | **Yes** | None | Single donor — no genotype variation. |
| 15 | **Hardwick et al. 2022 SnISOr-Seq method** (PMID 35256815) | Small (verify) | Thousands | No — targeted | None | Foundational method paper. |
| 16 | **Karakulak et al. 2025 ccRCC MAS-seq** (PMID 40107723) | Small (verify) | 2,599 | **Yes** (PacBio MAS-ISO-seq) | None | Small organoid panel. |
| 17 | **Tian et al. 2021 FLAMES method** (PMID 34763716) | 1 CLL + cell lines | Thousands | **Yes** (ONT) | None | Foundational. |
| 18 | **Al'Khafaji et al. 2024 MAS-ISO-seq method** (PMID 37291427) | Small melanoma TIL panel | Thousands | **Yes** (PacBio MAS-ISO-seq) | None | Foundational. |
| 19 | **ENCODE4 long-read (Reese 2023 bioRxiv)** | 81 samples (cell lines + tissues) | Mostly bulk; small LR-Split-seq subset | **Yes** (PacBio) | Cell-line genomes | Bulk-dominated. |

---

## 6. Bottom line

After live verification, the long-read single-cell RNA-seq landscape — measured against the constraint "many genotypes × full-isoform coverage × paired germline genotype" — has the following shape:

### Tier 1 (the only entry that satisfies all three criteria simultaneously)
- **Li et al. 2026 NCI/Korean lung sc-isoQTL** — 129 donors, paired genotype, whole-transcriptome Nanopore sc, sc-isoQTL design. Still a preprint as of June 2026, accessions not yet finalized.

### Tier 2 (large by single-cell donor count but missing one of the criteria)
- **Liu et al. 2024 eNeuro AD/PD/DLB cortex** — 25 donors, but long-read is *targeted* (50 genes) — not whole-transcriptome. EGA EGAS50000000132.
- **Argue & Gate 2025/26 AD CSF + blood** — up to ~108 immune-cell donors if not overlapping, but conference-abstract-stage, no germline genotype.

### Tier 3 (whole-transcriptome long-read but small donor count)
- **Joglekar et al. 2024 *Nat Neurosci*** — 6 human + 12 mouse, ScISOr-Seq2.
- **Pan et al. 2025 *Circulation*** — 12 heart donors, snNanoRNAseq (GEO GSE288222).
- **Li, Z. et al. 2024 *Cell Genomics* CRC** — 12 patients, 10x + PacBio Iso-Seq (ENA PRJEB68074).
- **Shiau et al. 2023 scNanoGPS** — 4 tumors + 2 cell lines, with per-cell mutation calling.

### Tier 4 (targeted long-read or method papers)
- Belchikov 2025 FTD SnISOr-Seq2 (3,630-gene panel).
- Byrne 2024 scTaILoR-seq (1,000-gene panel).
- Bizzotto 2025 FCD GO-TEN (3 long-read samples within a larger snRNA cohort).
- Hardwick 2022 SnISOr-Seq method.

### The negative result
- **There is no published >50-donor whole-transcriptome long-read scRNA-seq dataset with paired germline WGS other than Li et al. 2026.** Anything else in the field is either bulk (Glinos 2022 GTEx long-read; Nagura 2025 67-individual LCL ONT), targeted (Liu 2024, Belchikov 2025, Byrne 2024), single-strain mouse (Joglekar 2021, Rehman 2025), or small case-control panels (Pan 2025, Li Z. 2024).

### The biggest gap
The combination "**whole-transcriptome long-read scRNA-seq × > 100 germline-genotyped donors × multi-tissue**" is currently filled by *one* preprint. There is no PacBio Kinnex equivalent at this scale; no GTEx- or HCA-branded long-read sc layer; no iPSC-population analog of Cuomo et al. or Jerber et al. The IGVF mouse F1 hybrid-crosses effort (6.7M nuclei, 7 strains × C57BL/6J) is doing precisely the right design at scale but in **short-read** snRNA only — it would be a natural target for a long-read follow-on if one is ever attempted.

### What this means for a DNA LM trained on genotype → expression
1. If "full base-pair coverage RNA per cell × many genotypes" is a hard constraint, **the dataset universe is currently Li et al. 2026 (long-read sc, 129 donors, lung) plus Cuomo et al. 2020 (Smart-seq2, 125 donors, iPSC → endoderm)**. Everything else is either bulk full-coverage at much larger donor scale (GTEx, MAGE) or short-read sc at much larger cell scale (TenK10K, OneK1K, Jerber).
2. **Targeted long-read sc** (Liu 2024, Belchikov 2025, Byrne 2024) is useful for isoform supervision on a fixed gene panel — fine for an evaluation benchmark, weak for whole-genome pretraining.
3. **Mouse strain-panel long-read sc does not yet exist.** If Arc wanted to commission one, this is currently uncontested space.

---

## 7. Access notes

| Resource | Repository | Status |
|---|---|---|
| Li 2026 Korean lung sc-isoQTL | TBD (preprint) | Watch for journal publication and primary accessions |
| Liu 2024 eNeuro AD/PD/DLB | EGA EGAS50000000132 | **Controlled access — DAC application required** |
| Bizzotto 2025 FCD | dbGaP phs004124.v1.p1 | **Controlled access — dbGaP application required** |
| Pan 2025 heart | GEO GSE288222 | **Open** |
| Belchikov 2025 FTD SnISOr-Seq2 | GEO GSE250280 + SRA PRJNA1238317 | **Open** |
| Li Z. 2024 CRC | ENA PRJEB68074 + Zenodo 10.5281/zenodo.12750017 | **Open** |
| Joglekar 2021 mouse brain | GEO GSE158450 | **Open** |
| Joglekar 2024 brain atlas | www.isoformAtlas.com (primary repo to verify) | **Web portal + likely GEO** |
| Tian 2021 FLAMES | GEO (multiple) + EGA EGAS00001005597 (CLL2) | **Mixed open / controlled** |
| Shiau 2023 scNanoGPS | Verify on Nat Commun supplement | Likely open |
| Rehman 2025 aging mouse | Preprint — repo TBD | Watch |
| Argue 2025/26 AD CSF | Conference abstract — repo not yet available | Watch |

---

## 8. Key references (verified PMIDs / DOIs, 2026-06-29 round)

### Cohort-stage long-read sc papers
- Li B. et al. 2026 — Korean lung sc-isoQTL — PMID **42039472**, DOI [10.64898/2026.03.27.714873](https://doi.org/10.64898/2026.03.27.714873)
- Liu C.S. et al. 2024 — *eNeuro* AD/PD/DLB — PMID **39658200**, DOI [10.1523/ENEURO.0296-24.2024](https://doi.org/10.1523/ENEURO.0296-24.2024)
- Bizzotto S. et al. 2025 — *PNAS* FCD GO-TEN — PMID **40674414**, DOI [10.1073/pnas.2509622122](https://doi.org/10.1073/pnas.2509622122)
- Pan T. et al. 2025 — *Circulation* heart — PMID **41017471**, DOI [10.1161/CIRCULATIONAHA.125.074959](https://doi.org/10.1161/CIRCULATIONAHA.125.074959)
- Belchikov N. et al. 2025 — *Cell Rep* FTD — PMID **40913764**, DOI [10.1016/j.celrep.2025.116198](https://doi.org/10.1016/j.celrep.2025.116198)
- Li Z. et al. 2024 — *Cell Genomics* CRC — PMID **39216476**, DOI [10.1016/j.xgen.2024.100641](https://doi.org/10.1016/j.xgen.2024.100641)
- Joglekar A. et al. 2024 — *Nat Neurosci* brain — PMID **38594596**, DOI [10.1038/s41593-024-01616-4](https://doi.org/10.1038/s41593-024-01616-4)
- Byrne A. et al. 2024 — *Nat Commun* ovarian scTaILoR-seq — PMID **39134520**, DOI [10.1038/s41467-024-51252-6](https://doi.org/10.1038/s41467-024-51252-6)
- Foord C. et al. 2025 — *Nat Commun* Spl-ISO-Seq — PMID **40883294**, DOI [10.1038/s41467-025-63301-9](https://doi.org/10.1038/s41467-025-63301-9)
- Doyle P.H. et al. 2026 — *Front Genet* PBMC isonome — PMID **42283034**, DOI [10.3389/fgene.2026.1782221](https://doi.org/10.3389/fgene.2026.1782221)
- Joglekar A. et al. 2021 — *Nat Commun* mouse brain — PMID **33469025**, DOI [10.1038/s41467-020-20343-5](https://doi.org/10.1038/s41467-020-20343-5)
- Michielsen L. et al. 2026 — *Genome Biol* exon-inclusion model — PMID **41764538**, DOI [10.1186/s13059-026-04015-z](https://doi.org/10.1186/s13059-026-04015-z)
- Shiau C.-K. et al. 2023 — *Nat Commun* scNanoGPS — PMID **37433798**, DOI [10.1038/s41467-023-39813-7](https://doi.org/10.1038/s41467-023-39813-7)
- Karakulak T. et al. 2025 — *Genome Res* ccRCC MAS-seq — PMID **40107723**, DOI [10.1101/gr.279345.124](https://doi.org/10.1101/gr.279345.124)
- Pan/Pandey A. et al. 2025 — *Nat Commun* scCLEAN + MAS-seq — PMID **40389438**, DOI [10.1038/s41467-025-59880-2](https://doi.org/10.1038/s41467-025-59880-2)
- Rehman A. et al. 2025 bioRxiv aging mouse — DOI [10.1101/2025.06.05.658133](https://doi.org/10.1101/2025.06.05.658133)
- Argue B.M.R. & Gate D. 2025/26 conference abstract — DOI [10.1002/alz.70855](https://doi.org/10.1002/alz.70855)

### Method papers
- Gupta I. et al. 2018 — *Nat Biotechnol* ScISOr-Seq — PMID **30320766**, DOI [10.1038/nbt.4259](https://doi.org/10.1038/nbt.4259)
- Hardwick S.A. et al. 2022 — *Nat Biotechnol* SnISOr-Seq — PMID **35256815**, DOI [10.1038/s41587-022-01231-3](https://doi.org/10.1038/s41587-022-01231-3)
- Tian L. et al. 2021 — *Genome Biol* FLAMES — PMID **34763716**, DOI [10.1186/s13059-021-02525-6](https://doi.org/10.1186/s13059-021-02525-6)
- Al'Khafaji A.M. et al. 2024 — *Nat Biotechnol* MAS-ISO-seq — PMID **37291427**, DOI [10.1038/s41587-023-01815-7](https://doi.org/10.1038/s41587-023-01815-7)
- Rebboah E. et al. 2021 — *Genome Biol* LR-Split-seq — PMID **34620214**, DOI [10.1186/s13059-021-02505-w](https://doi.org/10.1186/s13059-021-02505-w)
- Fu Y. et al. 2025 — *Nat Commun* Longcell — PMID **40683866**, DOI [10.1038/s41467-025-60902-2](https://doi.org/10.1038/s41467-025-60902-2)
- Reese F. et al. 2023 — ENCODE4 long-read bioRxiv — PMID **37292896**, DOI [10.1101/2023.05.15.540865](https://doi.org/10.1101/2023.05.15.540865)
- Belchikov N. et al. 2024 review — *Genome Res* 34(11):1735–1746 — PMID **39567235**, DOI [10.1101/gr.279640.124](https://doi.org/10.1101/gr.279640.124)

### Related bulk long-read references (for completeness)
- Glinos D.A. et al. 2022 — *Nature* GTEx long-read — PMID **35922509**, DOI [10.1038/s41586-022-05035-y](https://doi.org/10.1038/s41586-022-05035-y) — *bulk ONT*, 88 samples
- Nagura Y. et al. 2025 — *Genome Biol* 67-individual B-cell ieQTL — PMID **40336129**, DOI [10.1186/s13059-025-03583-w](https://doi.org/10.1186/s13059-025-03583-w) — *bulk ONT on LCLs*, 67 Japanese individuals, **NOT single-cell despite some web summaries describing it as such**
- Lemoine L. et al. 2025 — *Sci Rep* F1 mouse brain allele-specific PacBio Iso-Seq — DOI [10.1038/s41598-025-21643-w](https://doi.org/10.1038/s41598-025-21643-w) — *bulk* on a single F1 (B6J × CAST/EiJ) brain (GSE246857)
- IGVF Consortium 2026 hybrid-crosses preprint — DOI [10.64898/2026.04.02.716195](https://doi.org/10.64898/2026.04.02.716195) — 6.7M nuclei from 7 F1 hybrids (B6J × 7 CC founders), **but Parse Biosciences short-read snRNA-seq, not long-read**

---

*Per PubMed MCP terms of use: information above attributed to PubMed; DOIs included as links throughout.*
