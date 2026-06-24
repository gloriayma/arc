# Jerber et al. 2021 — Data Access Log

**Paper:** Jerber*, Seaton*, Cuomo* et al., "Population-scale single-cell RNA-seq profiling across dopaminergic neuron differentiation," *Nature Genetics* 53:304–312, 2021.  
**DOI:** [10.1038/s41588-021-00801-6](https://doi.org/10.1038/s41588-021-00801-6)  
**PMID:** 33664506, **PMC:** PMC7610897  

---

## What this dataset is

- 215 human iPSC lines (HipSci consortium donors) differentiated toward midbrain dopaminergic neurons
- Profiled at 3 timepoints: **D11**, **D30**, **D52**
- >1 million cells total by scRNA-seq
- sc-eQTL mapped at each timepoint; 1,284 eQTL colocalize with neurological GWAS loci
- Paired WGS for all 215 donors (HipSci)

---

## ⚠️ Protocol note: this is NOT Smart-seq2

The `4_smartseq_datasets_with_genotype_variation.md` planning file lists this as Smart-seq2. **This is incorrect.** The GitHub repo (`single-cell-genetics/singlecell_neuroseq_paper`) explicitly contains a "10x Analysis Pipeline" built around CellRanger output. The ENA runs show paired-end Illumina HiSeq 4000 RNA-Seq, consistent with 10x Chromium droplet sequencing. The data is **3′-biased 10x Chromium, not full-length Smart-seq2**.

This matters for the DNA LM project: 10x data does NOT give full transcript coverage and does not capture splice junctions or UTRs. If full-length coverage is required, this dataset is not the right one (see Cuomo et al. E-MTAB-6945, which is also iPSC, 125 donors, Smart-seq2).

---

## Finding the correct accession (what didn't work)

### Attempt 1: Use E-MTAB-9610 from planning doc

The planning MD lists **E-MTAB-9610** as the accession. This is wrong.

```bash
curl -s "https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-9610"
```

Returns a completely unrelated study: "RNA-seq of co-culture experiment between human HUVEC endothelial cells and mouse glioma GL126 cells." Six samples, mouse only.

### Attempt 2: BioStudies text search

```bash
curl -s "https://www.ebi.ac.uk/biostudies/api/v1/search?query=Jerber+dopaminergic&pageSize=10"
```

Top hit: **S-EPMC7610897** — this is the PMC supplementary materials entry, not the primary data accession. It has supplementary figure JPEGs and an Excel table, not the sequencing data.

### Attempt 3: PubMed lookup

Used PubMed MCP tool to search for PMID. Found **PMID 33664506**. Then fetched the full PMC article text at `https://pmc.ncbi.nlm.nih.gov/articles/PMC7610897/` and extracted the Data Availability section. This worked and revealed the real accessions (see below).

---

## Correct accessions

| Repository | Accession | Access | Contents |
|---|---|---|---|
| **ENA** (open) | PRJEB38269 / ERP121676 | Open, no registration | 532 raw FASTQ runs (paired-end HiSeq 4000) |
| **EGA** (controlled) | Study EGAS00001002885, Dataset EGAD00001006157 | Managed access — requires application | Raw scRNA-seq reads |
| **Zenodo** (open) | [10.5281/zenodo.4333872](https://doi.org/10.5281/zenodo.4333872) | Open, CC-BY 4.0 | Processed count matrices (.h5 AnnData) + eQTL stats |

**Use the ENA or Zenodo routes** — no application needed for either.

---

## What's available without access application

### Zenodo (processed data) — start here

Direct download links (no login required):

| File | Size | Content |
|---|---|---|
| `D11.h5` | 5.88 GB | AnnData, day 11 cells |
| `D30.h5` | 5.68 GB | AnnData, day 30 cells |
| `D52.h5` | 13.05 GB | AnnData, day 52 cells |
| `all_timepoints_subsampled.h5` | 4.90 GB | 20% subsample across all timepoints |
| `eqtl_summary_stats.tar.gz` | 9.65 GB | eQTL summary statistics for 14 cell contexts |
| `coloc_neuroseq_25_traits.tsv.gz` | 9.88 MB | Colocalization: eQTL vs neurological GWAS |
| `coloc_gtex_25_traits.tsv.gz` | 62.77 MB | Colocalization: eQTL vs GTEx |

Download any file:
```bash
wget "https://zenodo.org/api/records/4333872/files/all_timepoints_subsampled.h5/content" \
     -O all_timepoints_subsampled.h5
```

Or for the smallest useful file first (the subsampled one at 4.9 GB):
```bash
wget -c "https://zenodo.org/api/records/4333872/files/all_timepoints_subsampled.h5/content" \
     -O all_timepoints_subsampled.h5
```

### ENA (raw FASTQ) — 532 runs, open access

```bash
# Get full list of runs with FTP paths
curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession%3DPRJEB38269&fields=run_accession,sample_accession,fastq_ftp&limit=532" > ena_runs.tsv

# Example single-run download (ERR4699948, ~9 GB paired):
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR469/008/ERR4699948/ERR4699948_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR469/008/ERR4699948/ERR4699948_2.fastq.gz
```

Each run is ~9–18 GB (paired FASTQ, ~140–290M read pairs). 532 runs total — this is the full raw data.

---

## Data format

The `.h5` files on Zenodo are **AnnData objects** saved in HDF5 format, loadable with scanpy:

```python
import scanpy as sc
adata = sc.read("D11.h5")
# adata.obs contains: cell barcode, celltype, time_point, donor line
# adata.X is the normalized count matrix
```

The notebooks in the GitHub repo confirm this:
- Repo: https://github.com/single-cell-genetics/singlecell_neuroseq_paper
- Load pattern: `adatafull = sc.read(file)` where file is a pre-processed per-timepoint or all-timepoints h5
- Cell type labels in `adata.obs`: FPP (floor plate progenitors), DA (dopaminergic neurons), NB (neuroblasts), Astro, Sert, Epen, etc.

---

## EGA controlled-access route (if raw reads are needed)

Raw reads also live at EGA under managed access:
- Study: EGAS00001002885
- Dataset: EGAD00001006157
- Application: https://ega-archive.org/studies/EGAS00001002885

This requires submitting a Data Access Request. For most purposes (expression analysis, not re-alignment), the Zenodo processed data is sufficient.

---

## Code and analysis resources

- GitHub repo: https://github.com/single-cell-genetics/singlecell_neuroseq_paper
  - `plotting_notebooks/` — Jupyter notebooks for all paper figures
  - `10x_analysis_pipeline/` — Snakemake pipeline (CellRanger → count matrix)
  - `differentiation_prediction_model/` — iPSC line differentiation outcome predictor
- Zenodo code archive: https://doi.org/10.5281/zenodo.4651413

---

## Recommended next step

Download `all_timepoints_subsampled.h5` (4.9 GB, 20% subsample across D11/D30/D52) for initial exploration:

```bash
cd /Users/gloria/dev/research/arc/sc_nona
mkdir -p data/jerber
wget -c "https://zenodo.org/api/records/4333872/files/all_timepoints_subsampled.h5/content" \
     -O data/jerber/all_timepoints_subsampled.h5
```

Then inspect with:
```python
import scanpy as sc
adata = sc.read("data/jerber/all_timepoints_subsampled.h5")
print(adata)
print(adata.obs.columns.tolist())
print(adata.obs['time_point'].value_counts())
```
