# Cuomo et al. 2020 — Data Access Log

**Paper:** Cuomo ASE, Seaton DD, McCarthy DJ et al., "Single-cell RNA-sequencing of differentiating iPS cells reveals dynamic genetic effects on gene expression," *Nature Communications* 11:810, 2020.  
**DOI:** [10.1038/s41467-020-14457-z](https://doi.org/10.1038/s41467-020-14457-z)  
**PMID:** 32041960, **PMC:** PMC7010688  

---

## What this dataset is

- 125 human iPSC lines (HipSci consortium) differentiated toward **definitive endoderm**
- 3 time points: day 0 (iPSC), day 1, day 3 of endoderm induction
- ~36,044 cells total after QC
- sc-eQTL mapped dynamically across differentiation stages
- Full-length Smart-seq2 — this is the real deal for the DNA LM project
- Paired genotype data via HipSci WGS

---

## ⚠️ Two errors in the planning doc (4_smartseq_datasets)

1. **Wrong accession**: E-MTAB-6945 points to an unrelated mouse neonate thymus T-cell study. Correct accessions: ENA PRJEB14362 (open), EGA EGAS00001002278 / EGAD00001005741 (controlled), Zenodo 3625024 (open).

2. **Wrong cell type**: The planning doc says "iPSC → cardiomyocyte differentiation." The paper is about **definitive endoderm** differentiation (iPSC → endoderm, 3-day protocol). Not cardiomyocytes.

---

## Finding the correct accession (what didn't work)

### Attempt 1: Use E-MTAB-6945 from planning doc

```bash
curl -s "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-6945"
```

Returns: "Single-cell analysis of T-cells from mouse neonate thymuses in a wild type control and a placental specific isoform of Igf2 knock-down (igf2-P0)." Completely unrelated.

### Attempt 2: PubMed search

Multiple PubMed searches with various keyword combinations failed to return any hits until using the author name "Cuomo Anna" + "quantitative trait loci" + "iPSC". That returned PMID 32041960.

### Attempt 3: PMC full text for data availability

Fetched `https://pmc.ncbi.nlm.nih.gov/articles/PMC7010688/` and extracted the Data Availability section. This revealed the correct accessions (see below).

### Attempt 4: ENA ERP016000 lookup

The paper lists ENA accession ERP016000 for the scRNA-seq data. Direct ENA count query returned 0 runs. However, querying ERP016000 as a secondary accession revealed it maps to primary BioProject **PRJEB14362**, which has 46,468 runs.

---

## Correct accessions

| Repository | Accession | Access | Contents |
|---|---|---|---|
| **ENA** (open) | PRJEB14362 (secondary: ERP016000) | Open, no registration | 46,468 raw FASTQ runs (Smart-seq2, HiSeq 2000, paired-end) |
| **EGA** (controlled) | Study EGAS00001002278, Dataset EGAD00001005741 | Managed access — requires application | Raw scRNA-seq reads |
| **Zenodo** (open) | [10.5281/zenodo.3625024](https://doi.org/10.5281/zenodo.3625024) | Open, CC-BY 4.0 | Processed count matrices |

---

## ⚠️ For base-pair resolution: use ENA raw FASTQs, NOT Zenodo

The Zenodo files contain **gene-level count matrices** — these are already summarized per gene per cell. You lose all sub-gene information: no splice junction reads, no UTR coverage, no allele-specific base calls, no positional signal. If the goal is BigWig-style full transcript coverage for a DNA language model, **Zenodo is the wrong thing to download**.

**What you need is the ENA raw FASTQs (PRJEB14362)**: these are the actual Smart-seq2 reads off the sequencer. Align them with STAR/HISAT2, then use bedtools/deeptools to generate per-base coverage BigWigs. This is the only route to base-pair resolution.

The pipeline: FASTQ → STAR alignment → BAM → bamCoverage (deeptools) → BigWig (base-pair coverage per cell).

---

## What's available without access application

### Zenodo (processed data) — gene counts only, NOT for base-pair resolution

| File | Size | Content |
|---|---|---|
| `raw_counts.csv.zip` | 2.54 GB | Raw gene × cell count matrix |
| `log_normalised_counts.csv.zip` | 2.48 GB | Log-normalized count matrix |
| `cell_metadata_cols.tsv` | 40.9 MB | Cell metadata (QC flags, time point, donor, etc.) |
| `ase_aggregated_by_donor_open_access_lines.tar.gz` | 844.7 MB | Allele-specific expression counts (HipSci open-access lines only) |

Download:
```bash
wget "https://zenodo.org/api/records/3625024/files/cell_metadata_cols.tsv/content" \
     -O cell_metadata_cols.tsv   # 40 MB, good first look

wget -c "https://zenodo.org/api/records/3625024/files/raw_counts.csv.zip/content" \
     -O raw_counts.csv.zip
```

**Format note:** The matrices are CSVs inside zip archives — not AnnData/h5ad. Rows = genes, columns = cells (or vice versa — check cell_metadata_cols.tsv first). You may want to load into AnnData after download:

```python
import pandas as pd
import anndata as ad
import numpy as np

counts = pd.read_csv("raw_counts.csv")     # ~36K cells × ~20K genes
meta   = pd.read_csv("cell_metadata_cols.tsv", sep="\t")
adata  = ad.AnnData(X=counts.values.T, obs=meta)
```

### ENA (raw FASTQs) — open access, 46K+ runs

```bash
# Get the full run manifest
curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession%3DPRJEB14362&fields=run_accession,sample_accession,fastq_ftp&limit=46468" > cuomo_ena_runs.tsv

# Example single-cell run (ERR1462954, ~1.8M paired reads, ~200 MB):
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR146/004/ERR1462954/ERR1462954_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR146/004/ERR1462954/ERR1462954_2.fastq.gz
```

Read counts per cell run: ~300K–2M read pairs, consistent with Smart-seq2 depth. Total raw data is ~46K runs × ~400 MB = very large; download selectively.

---

## Protocol confirmation

From the PMC methods section:

> "Single-cell transcriptomes of sorted cells were assayed as follows: reverse transcription and cDNA amplification was performed according to the SmartSeq2 protocol, and library preparation was performed using an Illumina Nextera kit."

This is **genuine full-length Smart-seq2**, run on Illumina HiSeq 2000. Full transcript coverage, informative for splice junctions, UTRs, allele-specific expression. This is the key advantage over Jerber (which is 10x Chromium).

---

## EGA controlled-access route (if raw reads needed)

- Study: EGAS00001002278
- Dataset: EGAD00001005741
- Application: https://ega-archive.org/studies/EGAS00001002278

HipSci data access also available via: http://www.hipsci.org

---

## Recommended next step

```bash
cd /Users/gloria/dev/research/arc/sc_nona
mkdir -p data/cuomo

# Start with metadata (40 MB, fast)
wget "https://zenodo.org/api/records/3625024/files/cell_metadata_cols.tsv/content" \
     -O data/cuomo/cell_metadata_cols.tsv

# Then get the raw counts (2.5 GB)
wget -c "https://zenodo.org/api/records/3625024/files/raw_counts.csv.zip/content" \
     -O data/cuomo/raw_counts.csv.zip
```
