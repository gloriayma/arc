# Smart-seq Datasets with Genotype Variation

*Compiled for Arc Institute scRNA-seq × DNA language model project.*
*Goal: datasets where multiple distinct individuals/genotypes are sequenced so a model can learn what DNA sequence variation does to transcriptional output.*

> **⚠ Major correction round — 2026-06-17.** The initial dataset list was generated from a Sonnet sub-agent using only training memory (no live PubMed/web access). Two rounds of verification via PubMed + ArrayExpress/GEO have found that *most entries were wrong* — incorrect protocols, accessions pointing to entirely unrelated studies, wrong journals, wrong authors, and in two cases the cited papers do not exist at all. **Net effect: the true Smart-seq2 + genotype-variation landscape is very sparse. Cuomo et al. 2020 is essentially the only large-scale (>50 donors) Smart-seq2 cohort with paired WGS that survives verification.**

### Correction log — round 1 (Cuomo, Jerber, Sarkar, Tung, Neavin, Krausgruber, Skelly, HCA, Tabula Sapiens, Villani)

> **Cuomo et al.:** E-MTAB-6945 was wrong. Correct accessions: ENA PRJEB14362 / ERP016000 (open), EGA EGAS00001002278 / EGAD00001005741 (controlled), Zenodo 10.5281/zenodo.3625024. Tissue: **definitive endoderm** (iPSC → endoderm, 3-day protocol), not cardiomyocyte. Protocol Smart-seq2 correct. Citation: *Nat Commun* 11:810 (2020), PMID 32041960. Notes: `doing_mds/2_cuomo_data.md`.

> **Jerber et al.:** (1) E-MTAB-9610 was wrong. Correct accessions: ENA PRJEB38269 (open), EGA EGAS00001002885 / EGAD00001006157 (controlled), Zenodo 10.5281/zenodo.4333872. (2) Protocol is **10x Chromium 3' v2**, not Smart-seq2. (3) Cells: **>1M**, not ~215K. Citation: *Nat Genet* 53:**304–312** (2021), PMID 33664506. **Coverage: 3'-biased, not full-length.** Notes: `doing_mds/1_jerber_data.md`.

> **Sarkar et al.:** Protocol is **Fluidigm C1 + UMIs (SMARTer chemistry)**, not Smart-seq2. Journal is ***PLoS Genet* 15(4):e1008045 (2019)**, not eLife. Accession is **GSE118723**, not GSE126030. Cells: 5,447. PMID 31002671.

> **Tung et al.:** Protocol is **Fluidigm C1 + UMIs**, not Smart-seq2. Citation: *Sci Rep* 7:**39921** (2017), PMID 28045081.

> **Neavin et al.:** Protocol is **10x Chromium**, not Smart-seq2. Journal ***Genome Biol* 22:76 (2021)**, not Cell Stem Cell. Donors: **110 total** (79 fibroblast + 31 iPSC). Cells: **~84K**. Accessions ArrayExpress **E-MTAB-9230** (fibroblast) / **E-MTAB-10060** (iPSC); GSE148780 was wrong. PMID 33673841. **3'-biased.**

> **Krausgruber et al.:** **DROPPED.** Mouse bulk multi-omics (epithelial/endothelial/fibroblast across 12 organs), not human scRNA-seq. ***Nature* 583:296 (2020)**. PMID 32612232.

> **Skelly et al. "mouse strain macrophages":** **DROPPED — paper does not exist** as described. GSE138852 was wrong (it's Grubman et al. AD snRNA-seq).

> **HCA pilot (Regev et al. 2017 *eLife*):** **DROPPED.** Perspective paper, not a dataset.

> **Tabula Sapiens:** Donor count **15** (v1), not 24 (24 = tissues).

> **Villani et al.:** Cells **2,422**.

### Correction log — round 2, 2026-06-17 (Hagemann-Jensen ×2, La Manno, Peng, Zhu/Zhong, Hodge, Kanton)

> **Hagemann-Jensen Smart-seq3 method paper:** Accession **E-MTAB-8735**, not GSE141630 (GSE141630 is an unrelated ovarian cancer study). Citation correct: *Nat Biotechnol* 38(6):708–714 (2020), PMID 32518404. Samples include mouse fibroblasts from a **CAST/EiJ × C57BL/6J F1 cross** (not single strain), HEK293FT, and a **5,376-cell HCA mixed-species benchmark** (not ~3,000 cells).

> **Hagemann-Jensen Smart-seq3xpress:** Wrong on multiple counts. Journal is ***Nat Biotechnol* 40(10):1452–1457 (2022)**, not *Genome Biol* 23:67. Accessions are **E-MTAB-11488 + E-MTAB-11452 + E-MTAB-11467**, not GSE173870 (GSE173870 is an unrelated EBF1 T-leukemia study). Cells: **26,260 hPBMCs** (not ~10,000). **PBMCs only, no fibroblasts.** Donors: 7. PMID 35637418.

> **La Manno et al. 2016 midbrain:** Protocol is **STRT-seq on Fluidigm C1 (UMI-based)**, NOT Smart-seq2. Citation/accession correct: *Cell* 167:566 (2016), GSE76381. 1,977 cells from 10 human embryos (6–11 GW). PMID 27716510.

> **Peng et al. pancreas 2019:** **DROPPED — fabricated citation.** No Peng et al. paper exists at *Cell Metab* 29:590 (2019). The real Peng et al. 2019 (***Cell Research* 29:725**, PMID 31273297) is a **PDAC tumor atlas** (10x Genomics, 57,530 cells across 24 PDAC + 11 control subjects) — not islets / T2D, not Smart-seq2, not 9 donors. GSE131886 also wrong (that accession is a different study, PMID 32354994).

> **"Zhu et al." cortex development — actually Zhong et al.:** First author is **Zhong** (Suijuan Zhong, Tang lab), not Zhu. Journal: ***Nature* 555(7697):524–528 (2018)***, not *Cell* 175:1575. Accession **GSE104276** is correct. Protocol Smart-seq2 confirmed. >2,300 cells from human PFC, GW 8–26 (no adult donors). PMID 29539641.

> **Hodge et al. cross-species cortical atlas (Allen):** Protocol is **SMART-Seq v4 Ultra Low Input on single nuclei (snRNA-seq)**, not Smart-seq2. Species: **Human + Mouse only — no macaque**. Cells: **15,928 nuclei from 8 human donors**, not ~76,000. Region: **middle temporal gyrus (MTG)**, not whole neocortex. Citation correct: *Nature* 573:61 (2019), PMID 31435019. Data on Allen Brain Map portal + NeMO Archive, no primary GEO.

> **Kanton et al. primate brain organoids:** Accessions are **ArrayExpress E-MTAB-7552 (10x) + E-MTAB-8234 (Fluidigm C1/Smart-seq2) + several others**, not GSE132672 (which is Bhaduri et al. cortical organoid stress, PMID 31996853). **Most data is 10x; only a smaller Fluidigm C1 subset is Smart-seq2.** Species: **Human + chimpanzee + macaque + bonobo** (no gorilla). Citation correct: *Nature* 574:418 (2019), PMID 31619793.

---

## 1. Verified Smart-seq2 (and Close Cousins) Cohorts with Genotype Variation

The set that survives both verification rounds. Smart-seq2 + multi-donor + paired WGS is essentially Cuomo et al. alone at meaningful scale.

| Name / Study | Citation | Accession | Protocol | Genotypes / Individuals | Cells | Tissue(s) | Species | Notes |
|---|---|---|---|---|---|---|---|---|
| **Cuomo et al. iPSC → endoderm eQTL** | *Nat Commun* 11:810 (2020), PMID 32041960 | ENA: PRJEB14362 (open); EGA: EGAD00001005741 (controlled); Zenodo: 10.5281/zenodo.3625024 (processed) | **Smart-seq2** | **125 donors** | ~36,044 | iPSC → definitive endoderm (day 0/1/3) | Human | **The de facto largest Smart-seq2 + genotype + paired-WGS cohort. Top priority.** |
| Picelli et al. Smart-seq2 method paper | *Nat Methods* 10:1096 (2013) | SRP028838 (SRA) | Smart-seq2 | ~6 donors (mixed) | ~100 cells | K562, fibroblasts, oocytes | Human | Foundational |
| Ziegenhain et al. protocol comparison | *Mol Cell* 65:631 (2017) | GSE75790 | Smart-seq2 + others | 3 cell lines (no genotype variation) | ~1,000 | Jurkat, K562, HEK293 | Human | Protocol benchmark |
| Tabula Muris (Smart-seq2 arm) | Schaum et al. 2018, *Nature* 562:367 (PMID 30283141) | GSE109774 | Smart-seq2 + 10x | **8 mice, all C57BL/6JN — no genotype variation** | 53,760 (SS2) + 70,118 (10x) | 20 organs | Mouse | Single strain |
| Tabula Muris Senis | Schaum et al. 2020, *Nature* 583:590 | GSE132040 | Smart-seq2 | ~18 mice, C57BL/6J — single strain | ~245,389 | 23 tissues | Mouse | Aging |
| **Tabula Sapiens (Smart-seq2 arm)** | Jones et al. 2022, *Science* 376:eabl4896 (PMID 35549404) | GSE201333 + CELLxGENE deposits | Smart-seq2 + 10x | **15 donors** (24 tissues; v1) | ~500K total (SS2 subset smaller) | 24 organs | Human | Multi-donor full-length component |
| Villani et al. DC atlas | *Science* 356:eaah4573 (2017), PMID 28428369 | GSE94820 | Smart-seq2 | Multiple healthy donors (exact count requires Table S1) | 2,422 | PBMCs (DCs, monocytes) | Human | High-resolution immune cell types |
| **"Zhong et al." (NOT Zhu) cortex development** | *Nature* 555(7697):524–528 (2018), PMID 29539641 | GSE104276 | **Smart-seq2** | Multiple fetal donors (GW 8–26; no adult) | >2,300 | Prefrontal cortex | Human | Developmental stages; fetal-only |
| Hagemann-Jensen Smart-seq3 method | *Nat Biotechnol* 38:708 (2020), PMID 32518404 | **E-MTAB-8735** | Smart-seq3 | Cell lines + CAST × B6 F1 mouse fibroblasts + HCA mixed-species | 5,376-cell HCA benchmark + 369 fibroblasts + HEK293FT | Multiple | Mixed | Method paper; CAST × B6 F1 has genotype variation but small scale |
| Hagemann-Jensen Smart-seq3xpress | *Nat Biotechnol* 40(10):1452–1457 (2022), PMID 35637418 | **E-MTAB-11488 + E-MTAB-11452 + E-MTAB-11467** | Smart-seq3xpress | **7 donors** | **26,260** | PBMCs only | Human | Cell-type isoform atlas; no paired WGS |
| Kanton et al. primate brain organoids (Smart-seq2/C1 subset only) | *Nature* 574:418 (2019), PMID 31619793 | E-MTAB-8234 (the SS2/C1 subset; full study spans many accessions) | **Mostly 10x; only smaller Fluidigm C1/SS2 subset** | 3–5 per species × **4 species (human, chimp, macaque, bonobo — NOT gorilla)** | ~750 organoid cells (full-length subset; total study is larger) | Brain organoids | Multi-species | Cross-species transcriptional divergence |

### Entries that were claimed as Smart-seq2 but use a different protocol

These are potentially useful for other purposes — just not Smart-seq2:

| Name / Study | Actual Protocol | Genotypes | Cells | Coverage | Notes |
|---|---|---|---|---|---|
| **Jerber et al. neurons** (ENA: PRJEB38269) | **10x Chromium 3' v2** | 215 donors | >1M | 3'-biased | Largest iPSC sc-eQTL overall |
| **Neavin et al.** (E-MTAB-9230, E-MTAB-10060) | **10x Chromium** | 110 donors total | ~84K | 3'-biased | *Genome Biol* 22:76 (2021) |
| **Sarkar et al. Yoruba iPSC** (GSE118723) | **Fluidigm C1 + UMIs** | 53 donors | 5,447 | Partial (5'/3' bias) | *PLoS Genet* 2019 |
| **Tung et al. Yoruba iPSC** (GSE77288) | **Fluidigm C1 + UMIs** | 3 lines × 3 reps | ~864 | Partial | *Sci Rep* 7:39921 (2017) |
| **La Manno et al. midbrain** (GSE76381) | **STRT-seq on Fluidigm C1 (UMI-based)** | 10 human embryos | 1,977 (human subset) | 5'-end UMI (not full-length) | *Cell* 167:566 (2016) |
| **Hodge et al. (Allen) cortical** | **SMART-Seq v4 Ultra Low Input snRNA-seq** (single nuclei) | 8 human donors | 15,928 nuclei | Full-length nucleus reads | *Nature* 573:61 (2019); human + mouse only, MTG region |

### Entries fully removed (fabricated or wrong study)

- **Peng et al. pancreas 2019 *Cell Metab*** — paper does not exist. Real Peng 2019 is *Cell Research* 29:725, a 10x PDAC tumor atlas.
- **Krausgruber et al. 2020** — mouse bulk multi-omics, not human scRNA.
- **Skelly et al. "mouse strain macrophages"** — paper does not exist as described.
- **HCA pilot (Regev et al. 2017)** — perspective paper, not a dataset.

---

## 2. Large Population-Scale scRNA-seq (Genotype × Expression at Scale — Not Full Coverage)

Not Smart-seq, but premier resources for genotype-linked transcriptional variation at population scale.

| Name / Study | Citation | Accession | Protocol | Donors | Cells | Tissue(s) | Species |
|---|---|---|---|---|---|---|---|
| **OneK1K** | Yazar et al. 2022, *Science* 376:eabf3041 | EGA: EGAS00001005634 | 10x Chromium v3 | 982 | ~1.27M | PBMCs | Human |
| **TenK10K Phase 1** (Garvan / Yazar lab) | Preprints: Cuomo et al. *medRxiv* 2025 (10.1101/2025.03.20.25324352); no journal pub yet | EGAS50000001654 (WGS) + EGAS50000001653 (sc) | 10x | 1,925 (WGS-matched) | ~5M | PBMCs | Human |
| Blueprint sc-eQTL | van der Wijst et al. 2020, *Nat Genet* 52:1194 | EGA: EGAD00001006114 | 10x | 120 | ~1.3M | PBMCs | Human |
| **GTEx snRNA-seq** | GTEx Consortium 2020, *Science* 369:1318, PMID 32913098 | dbGaP phs000424 | 10x snRNA-seq | ~838 (analyzed) | ~500K nuclei | 49 tissues | Human |
| iPSCORE | **Panopoulos et al. 2017, *Stem Cell Reports* 8:1086–1100**, PMID 28410642 (the prior *Cell Stem Cell* citation was wrong) | dbGaP phs001325 | Bulk + sc subset | 222 | Bulk primarily | iPSC | Human |
| CZI CELLxGENE Census | CZI 2023 | cellxgene.cziscience.com | Mixed | >1,000 aggregated | >50M | All tissues | Human, Mouse |
| eQTLGen | Võsa et al. 2021, *Nat Genet* 53:1300, PMID 34475573 | eqtlgen.org | Bulk RNA-seq meta | 31,684 (meta-analysis) | N/A (bulk) | Whole blood | Human |

For the full-coverage × genotype landscape beyond Smart-seq, see `4.5_biggest_full_coverage_scrna_datasets_with_genotype_variation.md`.

---

## 3. Multi-Species Datasets

| Name / Study | Citation | Accession | Protocol | Species | Notes |
|---|---|---|---|---|---|
| Kanton et al. primate brain organoids | *Nature* 574:418 (2019), PMID 31619793 | E-MTAB-7552 (10x), E-MTAB-8234 (C1/SS2), + more | Mostly 10x; small SS2/C1 subset | Human, chimp, macaque, **bonobo** | Cross-species; **no gorilla** |
| Hodge et al. cortical (Allen) | *Nature* 573:61 (2019), PMID 31435019 | Allen Brain Map / NeMO | SMART-Seq v4 snRNA-seq | **Human + Mouse** (no macaque) | MTG region; ~16K nuclei |
| La Manno et al. midbrain | *Cell* 167:566 (2016), PMID 27716510 | GSE76381 | STRT-seq on Fluidigm C1 | Human, Mouse | Developmental; 10 human embryos |
| Cardoso-Moreira et al. evo transcriptome | *Nature* 571:505 (2019) | ERP109255 | Bulk RNA-seq | 7 species | Bulk only |

---

## 4. Mouse Strain Panels

No verified Smart-seq2 mouse-strain single-cell panel. Real strain resources are bulk:

| Name / Study | Protocol | Strains | Notes |
|---|---|---|---|
| Collaborative Cross founder bulk/sc | Mixed | 8 CC founders | Paired strain genomes |
| Diversity Outbred (Jackson Lab) | Mostly bulk; some 10x | >100 DO mice | Paired WGS reconstruction |
| BXD panel | Bulk RNA-seq | >140 BXD strains | GeneNetwork.org |
| Sanger Mouse Genomes Project + scRNA-seq | 10x mostly | 17 inbred strains | Full WGS on all 17 strains |
| ImmGen | Bulk microarray/RNA-seq | >200 mouse strains | Foundational |

---

## 5. Revised Priority Ranking

For training a DNA LM on genotype → expression with full-coverage signal.

| Rank | Dataset | Protocol | Genotypes | Cells | Paired WGS? | Notes |
|---|---|---|---|---|---|---|
| 1 | **Cuomo et al. iPSC → endoderm** (PRJEB14362) | **Smart-seq2** | 125 donors | ~36K | Yes | The clear top Smart-seq2 + genotype + WGS dataset. Full-length coverage. |
| 2 | **Jerber et al. neurons** (PRJEB38269) | **10x 3' v2** | 215 donors | >1M | Yes | Largest iPSC sc-eQTL cohort, but **3'-biased — not full coverage** |
| 3 | **TenK10K Phase 1** (EGAS50000001654 + EGAS50000001653) | 10x | 1,925 | ~5M | Yes | Largest sc-eQTL by donor count; 3'-biased |
| 4 | **OneK1K** (EGAS00001005634) | 10x | 982 | 1.27M | Yes | Largest sc-eQTL by cell count; 3'-biased |
| 5 | **GTEx snRNA-seq** (phs000424) | 10x snRNA | ~838 | ~500K nuclei | Yes | Multi-tissue; 3'-biased |
| 6 | **Tabula Sapiens (SS2 arm)** (GSE201333) | Smart-seq2 + 10x | 15 donors | ~500K total | Partial | Multi-organ Smart-seq2 arm |
| 7 | **Sarkar et al. Yoruba iPSC** (GSE118723) | Fluidigm C1 + UMIs | 53 donors | 5,447 | Yes | Partial coverage |
| 8 | **Neavin et al.** (E-MTAB-9230, E-MTAB-10060) | 10x | 110 donors | ~84K | Partial | 3'-biased |
| 9 | **Hagemann-Jensen Smart-seq3xpress** (E-MTAB-11488 et al.) | Smart-seq3xpress | 7 donors | 26,260 | No | Full-length + UMI; no WGS |
| 10 | **Zhong et al. fetal PFC** (GSE104276) | Smart-seq2 | Multiple GW8–26 fetal | >2,300 | No | Developmental |
| 11 | **Kanton et al. primate organoids** (E-MTAB-7552 + 8234) | Mostly 10x; SS2/C1 subset | 3–5 per species × 4 | ~750 (SS2 subset) | Genomes available | Cross-species (human/chimp/macaque/bonobo) |
| 12 | **Tung et al. Yoruba iPSC** (GSE77288) | Fluidigm C1 + UMIs | 3 × 3 reps | ~864 | No (imputed) | Population genetics design |

---

## 6. Bottom Line

After two verification rounds, the Smart-seq2 + genotype-variation landscape that survives is **essentially one cohort: Cuomo et al. 2020** (125 donors, iPSC → endoderm, ~36K cells, paired WGS, ENA PRJEB14362).

For DNA LM training:

1. **Cuomo et al. 2020** — only large Smart-seq2 + genotype + WGS cohort. Use for full-coverage signal at single-cell resolution.
2. **Jerber et al. 2021** — much larger (215 donors, >1M cells, paired WGS) but **10x 3'-biased**.
3. **TenK10K Phase 1** — largest by donor count (1,925) but 10x, no peer review yet.
4. **Bulk full-coverage cohorts** (GTEx, MAGE, Geuvadis, etc.) — see file 4.5 — dominate single-cell when single-cell resolution is not required.
5. **Long-read single-cell** (Li et al. 2026 Korean lung, Nanopore, 129 donors) — see file 4.5 — first population-scale long-read sc-isoQTL.

---

## 7. Access Notes

Data access applications likely needed for:
- **GTEx** (phs000424): dbGaP
- **OneK1K / TenK10K**: EGA
- **Blueprint sc-eQTL**: EGA
- **Cuomo et al.** controlled tier: EGA EGAS00001002278 (processed counts open on Zenodo 10.5281/zenodo.3625024)
- **Jerber et al.** controlled tier: EGA EGAS00001002885 (processed counts open on Zenodo 10.5281/zenodo.4333872)

### Primary verification sources (2026-06-17 rounds)

PubMed PMIDs / DOIs used for verification:
- Cuomo et al. 2020 — PMID 32041960, DOI 10.1038/s41467-020-14457-z
- Jerber et al. 2021 — PMID 33664506, DOI 10.1038/s41588-021-00801-6
- Sarkar et al. 2019 — PMID 31002671, DOI 10.1371/journal.pgen.1008045
- Tung et al. 2017 — PMID 28045081, DOI 10.1038/srep39921
- Neavin et al. 2021 — PMID 33673841, DOI 10.1186/s13059-021-02293-3
- Krausgruber et al. 2020 — PMID 32612232, DOI 10.1038/s41586-020-2424-4
- Tabula Sapiens 2022 — PMID 35549404, DOI 10.1126/science.abl4896
- Tabula Muris 2018 — PMID 30283141, DOI 10.1038/s41586-018-0590-4
- Villani et al. 2017 — PMID 28428369, DOI 10.1126/science.aah4573
- Hagemann-Jensen Smart-seq3 2020 — PMID 32518404, DOI 10.1038/s41587-020-0497-0
- Hagemann-Jensen Smart-seq3xpress 2022 — PMID 35637418, DOI 10.1038/s41587-022-01311-4
- La Manno 2016 — PMID 27716510, DOI 10.1016/j.cell.2016.09.027
- Peng 2019 (real, not the fabricated one) — PMID 31273297, DOI 10.1038/s41422-019-0195-y
- Zhong 2018 — PMID 29539641, DOI 10.1038/nature25980
- Hodge 2019 — PMID 31435019, DOI 10.1038/s41586-019-1506-7
- Kanton 2019 — PMID 31619793, DOI 10.1038/s41586-019-1654-9
- Panopoulos 2017 (iPSCORE, correct citation) — PMID 28410642, DOI 10.1016/j.stemcr.2017.03.012
