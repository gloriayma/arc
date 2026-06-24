# Data Download Guide — Cuomo et al. 2020 (Cuomo scRNA + HipSci genomes)

**Goal:** Build per-donor N-masked personalized reference genomes, align Smart-seq2 reads with WASP bias correction, split into haplotype-resolved BAMs with SNPsplit, and produce per-cell per-haplotype base-pair coverage BigWigs for DNA LM training.

**Access status:** Everything here is open access. No EGA application required.  
**Note on 21 missing donors:** donors cicb, dixh, eoxi, fasu, fejf, guss, guyj, hegp, koqx, nocf, nudd, oebj, oibg, oojs, pulk, tavh, tout, ueah, veku, walu, zoio are managed access (EGA only). See `doing_mds/2.2b_access_verification.md` for details.

> **Pipeline change from earlier draft (2026-06-22):** The previous version built full diploid FASTAs (two haplotype sequences per donor) and split haplotypes by chromosome name after alignment. This was wrong for two reasons: (1) aligning short reads to a diploid concatenated reference causes multi-mapping explosion — reads with no overlapping variant map equally to both haplotypes; (2) the coordinate renaming (`sed 's/_h2\b//g'`) breaks if any indels were applied, since indels shift positions downstream. The correct approach (from `doing_mds/2.4_diploid_alignment.md`) is: **N-masked single reference + STAR+WASP + SNPsplit**. This avoids both problems. See 2.4 for full rationale.

---

## Overview of what to download

| What | Size | Why |
|---|---|---|
| Cell metadata (Zenodo) | 41 MB | Maps cell IDs → donor names; needed to group cells by donor |
| ENA CRAM manifest | <1 MB | Maps all 46,468 run accessions → CRAM download URLs |
| HipSci WGS manifest | <1 MB | Maps donor names → WGS VCF download URLs |
| hs37d5 reference genome | ~3 GB | Required to decode CRAMs + base for personalized references |
| GENCODE v19 GTF | ~35 MB | Gene annotation for STAR splice junction index (one file, shared across all donors) |
| HipSci imputed BCFs (genotypes) | ~6 GB | Per-donor phased genotypes for N-masking and SNPsplit |
| scRNA-seq CRAMs (~31 GB open-access subset) | ~31 GB | The actual single-cell RNA-seq reads |

**Total download: ~40 GB**

---

## Step 0 — Create a project directory

```bash
mkdir -p cuomo_data/{metadata,metadata/cells_by_donor,metadata/cram_urls_by_donor,\
references,references/personalized,genotypes,genotypes/imputed_bcf,genotypes/per_donor,\
crams,bams,bigwigs}
cd cuomo_data
```

---

## Step 1 — Download metadata (start here, fast)

These small files tell you which cells belong to which donors, which is needed to group cells by donor for alignment.

### 1a. Cell metadata from Zenodo

```bash
wget -c "https://zenodo.org/api/records/3625024/files/cell_metadata_cols.tsv/content" \
     -O metadata/cell_metadata_cols.tsv
# Size: 41 MB. Has columns: cell_id, donor, time_point, experiment, plate, etc.
```

### 1b. ENA CRAM manifest (all 46,468 runs)

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB14362&result=read_run&fields=run_accession,sample_accession,library_name,submitted_ftp&format=tsv&limit=50000" \
     > metadata/ena_cram_manifest.tsv
# Returns 46,468 rows. submitted_ftp column has the CRAM + CRAI URLs.
```

### 1c. HipSci WGS manifest

```bash
curl -s "https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/archive_datasets/ENA.ERP017015.wgs.normals.analysis_files.tsv" \
     > metadata/hipsci_wgs_manifest.tsv
# Maps HipSci line names → WGS gVCF download URLs (530 samples)
```

### 1d. Identify open-access cells and group by donor

```python
import pandas as pd

meta = pd.read_csv("metadata/cell_metadata_cols.tsv", sep="\t")

# 21 managed-access donor short names
managed = {'cicb','dixh','eoxi','fasu','fejf','guss','guyj','hegp',
           'koqx','nocf','nudd','oebj','oibg','oojs','pulk','tavh',
           'tout','ueah','veku','walu','zoio'}

# Donor column contains e.g. "HPSI0114i-bezi_1" — extract the short name
meta['donor_short'] = meta['donor'].str.extract(r'HPSI\d+i-([a-z]+)_?\d*')[0]
open_meta = meta[~meta['donor_short'].isin(managed)]

print(f"Open-access cells: {len(open_meta)} / {len(meta)}")

# Save per-donor cell lists (needed to group CRAMs by donor)
import os
os.makedirs("metadata/cells_by_donor", exist_ok=True)
for donor, group in open_meta.groupby('donor'):
    group.to_csv(f"metadata/cells_by_donor/{donor}.tsv", sep="\t", index=False)

open_meta.to_csv("metadata/open_access_cells.tsv", sep="\t", index=False)
```

---

## Step 2 — Download reference files

### 2a. hs37d5 reference genome (REQUIRED to decode CRAMs, and base for personalized refs)

The CRAMs are compressed against this exact reference. You must have it to decode them. It is also the base sequence that per-donor variants are applied on top of to build personalized references.

```bash
wget -c "https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/reference/hs37d5.fa.gz" \
     -O references/hs37d5.fa.gz
# Size: ~3 GB. GRCh37 + decoy sequences + EBV + mitochondrial rCRS.

gunzip references/hs37d5.fa.gz
samtools faidx references/hs37d5.fa
```

### 2b. GENCODE v19 gene annotation (GRCh37, matches the paper)

One GTF, shared across all donors — no duplication needed.

```bash
wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz" \
     -O references/gencode.v19.annotation.gtf.gz
# Size: ~35 MB
gunzip references/gencode.v19.annotation.gtf.gz
```

---

## Step 3 — Download genotypes and build N-masked personalized FASTAs

For each donor, we build a single personalized reference FASTA where:
- Positions where the donor's two copies agree and differ from the standard reference: substituted with the donor's allele
- Positions where the donor's two copies disagree (heterozygous): replaced with `N`

The N-masking means the aligner cannot favour either allele at heterozygous sites. WASP (Step 5b) then catches any residual bias from reads spanning multiple het sites.

> **Why not build two full haplotype FASTAs?** Short reads (~150 bp) that don't overlap any variant site map equally to both haplotypes, producing duplicate alignments throughout the genome. STAR's junction database also doesn't transfer cleanly across haplotypes. The N-masked single-reference approach avoids both problems and is the field standard for short-read RNA-seq. Haplotype splitting is handled after alignment by SNPsplit (Step 5c), which uses the phasing information from the VCF to assign reads that do overlap variant sites. See `doing_mds/2.4_diploid_alignment.md` §2 and §3 for benchmarks.

> **What about indels?** The N-masking approach only handles SNPs cleanly. We deliberately exclude indels: applying a donor-specific insertion or deletion shifts all downstream coordinates, breaking compatibility with the GTF. For this project, indels in expressed regions are treated as a known limitation — the aligner will soft-clip reads spanning donor-specific indels rather than placing them perfectly. This is consistent with current best practice in the field.

### 3a. Download imputed joint BCF (all open-access donors, ~6 GB)

The imputed BCF is statistically phased (`0|1` format, SHAPEIT v2), which is what allows SNPsplit to later assign reads to hap1 vs hap2. It covers ~6–10M common variants per donor.

```bash
for CHR in {1..22} X; do
    wget -c "https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/data/vep_openaccess_bcf/chr${CHR}.bcf" \
         -O genotypes/imputed_bcf/chr${CHR}.bcf
    wget -c "https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/data/vep_openaccess_bcf/chr${CHR}.bcf.csi" \
         -O genotypes/imputed_bcf/chr${CHR}.bcf.csi
done
# Total: ~6 GB
```

### 3b. Extract per-donor VCF

```bash
DONOR="HPSI0114i-bezi_1"   # replace with actual donor name

# Extract this donor from each chromosome and concatenate
for CHR in {1..22} X; do
    bcftools view -s ${DONOR} genotypes/imputed_bcf/chr${CHR}.bcf \
        -Oz -o genotypes/per_donor/${DONOR}.chr${CHR}.vcf.gz
    bcftools index -t genotypes/per_donor/${DONOR}.chr${CHR}.vcf.gz
done

bcftools concat \
    $(for CHR in {1..22} X; do echo genotypes/per_donor/${DONOR}.chr${CHR}.vcf.gz; done) \
    -Oz -o genotypes/per_donor/${DONOR}.vcf.gz
bcftools index -t genotypes/per_donor/${DONOR}.vcf.gz

rm genotypes/per_donor/${DONOR}.chr*.vcf.gz*
```

### 3c. Build N-masked personalized FASTA (replaces the diploid FASTA approach)

```bash
DONOR="HPSI0114i-bezi_1"
VCF="genotypes/per_donor/${DONOR}.vcf.gz"

# ── Step 1: Apply homozygous-alt SNPs (personalize the reference) ──────
# These are positions where both copies agree and differ from hs37d5.
# We substitute the donor's allele directly — no ambiguity.
bcftools view -g hom -e 'GT="0/0" || GT="0|0"' -V indels ${VCF} \
    -Oz -o genotypes/per_donor/${DONOR}.hom_alt_snps.vcf.gz
bcftools index -t genotypes/per_donor/${DONOR}.hom_alt_snps.vcf.gz

bcftools consensus \
    -f references/hs37d5.fa \
    genotypes/per_donor/${DONOR}.hom_alt_snps.vcf.gz \
    > references/personalized/${DONOR}.hom_applied.fa

# ── Step 2: Extract het SNP positions and N-mask them ──────────────────
# These are positions where the two copies differ. We write N so the aligner
# treats both alleles equally.
bcftools view -g het -V indels ${VCF} \
    | bcftools query -f '%CHROM\t%POS0\t%POS\n' \
    > references/personalized/${DONOR}.het_snps.bed

bedtools maskfasta \
    -fi references/personalized/${DONOR}.hom_applied.fa \
    -bed references/personalized/${DONOR}.het_snps.bed \
    -fo references/personalized/${DONOR}.Nmasked.fa

samtools faidx references/personalized/${DONOR}.Nmasked.fa

# ── Step 3: Prepare SNPsplit input file (do once per donor) ───────────
# SNPsplit needs a file listing het SNP positions with their phasing.
# Format: SNP-ID  CHROM  POS  STRAND  REF/ALT
bcftools query -s ${DONOR} \
    -i 'GT[*]="het" && TYPE="snp"' \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' \
    ${VCF} \
    > genotypes/per_donor/${DONOR}.snpsplit_input.txt
# SNPsplit reads the phase (0|1 vs 1|0) from the GT column to know which
# allele is hap1 vs hap2.

# Clean up intermediate
rm references/personalized/${DONOR}.hom_applied.fa
```

### 3d. Build donor STAR index (~30 GB, one standard GTF)

```bash
STAR_IDX="references/STAR_personalized/${DONOR}"
mkdir -p ${STAR_IDX}

STAR --runMode genomeGenerate \
     --genomeDir ${STAR_IDX} \
     --genomeFastaFiles references/personalized/${DONOR}.Nmasked.fa \
     --sjdbGTFfile references/gencode.v19.annotation.gtf \
     --sjdbOverhang 74 \
     --runThreadN 16
# ~30 min, ~30 GB disk (one standard GTF — no duplication for hap2 needed)
```

---

## Step 4 — Download the scRNA-seq CRAMs

Each CRAM is one cell. There are 46,468 total; we want only the open-access donors' cells.

### 4a. Generate per-donor CRAM URL lists

```python
import pandas as pd

ena = pd.read_csv("metadata/ena_cram_manifest.tsv", sep="\t")
open_cells = pd.read_csv("metadata/open_access_cells.tsv", sep="\t")

# Inspect the key linking ENA runs to cell metadata
print(ena[['run_accession','library_name']].head())
print(open_cells[['cell_id']].head())

# Join on library_name / cell_id (adjust column names after inspection)
open_ena = ena[ena['library_name'].isin(open_cells['cell_id'])]
open_ena = open_ena.merge(open_cells[['cell_id','donor']], left_on='library_name', right_on='cell_id')

import os
os.makedirs("metadata/cram_urls_by_donor", exist_ok=True)
for donor, group in open_ena.groupby('donor'):
    urls = group['submitted_ftp'].str.split(';').explode()
    urls = urls[urls.str.endswith('.cram')]
    urls.to_csv(f"metadata/cram_urls_by_donor/{donor}.txt", index=False, header=False)
```

### 4b. Download CRAMs for one donor

```bash
DONOR="HPSI0114i-bezi_1"
mkdir -p crams/${DONOR}

cat metadata/cram_urls_by_donor/${DONOR}.txt \
    | xargs -P 8 -I{} wget -c -q {} -P crams/${DONOR}/

# Download CRAI index files
sed 's/\.cram$/.cram.crai/' metadata/cram_urls_by_donor/${DONOR}.txt \
    | xargs -P 8 -I{} wget -c -q {} -P crams/${DONOR}/
```

---

## Step 5 — Pipeline: CRAM → FASTQ → STAR+WASP → SNPsplit → BigWig

Run this for each cell, after building the donor's STAR index (Step 3d).

```bash
DONOR="HPSI0114i-bezi_1"
CELL="19776_1#2"              # replace with actual cell ID
CRAM="crams/${DONOR}/${CELL}.cram"
STAR_IDX="references/STAR_personalized/${DONOR}"
REF="references/hs37d5.fa"   # needed to decode CRAM

mkdir -p bams/${DONOR}/${CELL} bigwigs/${DONOR}

# ── Step 5a: CRAM → FASTQ ──────────────────────────────────────────────
# Discards the BWA-MEM alignment, recovers raw reads losslessly.
samtools collate -uOn 128 --reference ${REF} ${CRAM} /tmp/${CELL} \
    | samtools fastq -F 0xB00 \
                     -1 /tmp/${CELL}_R1.fq.gz \
                     -2 /tmp/${CELL}_R2.fq.gz \
                     -0 /dev/null -s /dev/null -

# ── Step 5b: STAR alignment with WASP bias correction ─────────────────
# --varVCFfile: tells STAR which positions are het in this donor
# --waspOutputMode SAMtag: tags each read vW:i:1 (pass) or vW:i:2+ (fail)
#   WASP logic: at each het site, flip the allele and re-align; if the read
#   maps differently, it was reference-biased → discard.
# --alignEndsType EndToEnd: required by SNPsplit (cannot handle soft-clipped reads)

STAR --runMode alignReads \
     --genomeDir ${STAR_IDX} \
     --readFilesIn /tmp/${CELL}_R1.fq.gz /tmp/${CELL}_R2.fq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI NM MD vA vG vW \
     --twopassMode Basic \
     --varVCFfile genotypes/per_donor/${DONOR}.vcf.gz \
     --waspOutputMode SAMtag \
     --alignEndsType EndToEnd \
     --runThreadN 4 \
     --outFileNamePrefix bams/${DONOR}/${CELL}/

BAM="bams/${DONOR}/${CELL}/Aligned.sortedByCoord.out.bam"
samtools index ${BAM}

# ── Step 5c: Filter WASP-passing reads, then split by haplotype ────────
# 1. Keep only WASP-passing reads (vW:i:1 tag)
samtools view -h ${BAM} \
    | awk '/^@/ || /vW:i:1/' \
    | samtools view -b - \
    > bams/${DONOR}/${CELL}/wasp_filtered.bam
samtools index bams/${DONOR}/${CELL}/wasp_filtered.bam

# 2. SNPsplit: assign reads to hap1 or hap2 using the phased genotype.
#    Reads overlapping a het SNP where the read carries the '0|1' allele
#    go to genome1 (hap1); reads carrying the '1|0' allele go to genome2 (hap2).
#    Reads with no overlapping het SNP (~60-80% of reads) go to unassigned.bam.
SNPsplit --paired \
         --snp_file genotypes/per_donor/${DONOR}.snpsplit_input.txt \
         bams/${DONOR}/${CELL}/wasp_filtered.bam
# Outputs (in same directory as input):
#   wasp_filtered.genome1.bam  → hap1-assigned reads
#   wasp_filtered.genome2.bam  → hap2-assigned reads
#   wasp_filtered.unassigned.bam → no overlapping het SNP (not phase-assignable)

samtools index bams/${DONOR}/${CELL}/wasp_filtered.genome1.bam
samtools index bams/${DONOR}/${CELL}/wasp_filtered.genome2.bam

# ── Step 5d: BAM → BigWig (one per haplotype + one diploid-summed) ─────
# hap1 track
bamCoverage \
    --bam bams/${DONOR}/${CELL}/wasp_filtered.genome1.bam \
    --outFileName bigwigs/${DONOR}/${CELL}.h1.bw \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 4

# hap2 track
bamCoverage \
    --bam bams/${DONOR}/${CELL}/wasp_filtered.genome2.bam \
    --outFileName bigwigs/${DONOR}/${CELL}.h2.bw \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 4

# Diploid-summed track (all WASP-passing reads regardless of phase assignment)
# Useful if you want a single expression track per cell without haplotype splitting.
bamCoverage \
    --bam bams/${DONOR}/${CELL}/wasp_filtered.bam \
    --outFileName bigwigs/${DONOR}/${CELL}.diploid.bw \
    --binSize 1 \
    --normalizeUsing CPM \
    --numberOfProcessors 4

# Clean up temporaries
rm /tmp/${CELL}_R1.fq.gz /tmp/${CELL}_R2.fq.gz
```

> **Caveat on haplotype tracks:** ~60–80% of Smart-seq2 reads will land in `unassigned.bam` because they don't overlap any het SNP. The h1 and h2 BigWigs therefore have sparser coverage than the diploid-summed track. For DNA LM training, the diploid-summed track may be more practical as a primary signal, with haplotype tracks available for fine-grained allele-specific analysis.

> **Training data pairing:** For the DNA LM, a training example is: donor's hs37d5-based personalized sequence window → per-cell coverage BigWig aligned to those same coordinates. The BigWigs are in hs37d5 coordinate space (matching the personalized reference and the VCF). Comparing coverage across donors requires mapping to a common coordinate system, but within a donor everything is internally consistent.

---

## Step 6 — Processing order (per donor)

Process one donor at a time to manage disk space:

```
1. Step 3b: extract per-donor VCF                   (~5 min)
2. Step 3c: build N-masked FASTA + SNPsplit input   (~10 min)
3. Step 3d: build STAR index                        (~30 min, ~30 GB disk)
4. Step 4b: download CRAMs                          (depends on # cells)
5. Step 5:  align + WASP filter + SNPsplit          (~3–7 min/cell × ~350 cells ≈ ~20 hrs on 4 cores)
6. Delete STAR index                                (free 30 GB before next donor)
```

---

## Summary of accessions

| Data | Accession | URL |
|---|---|---|
| scRNA-seq CRAMs | ENA PRJEB14362 | https://www.ebi.ac.uk/ena/browser/view/PRJEB14362 |
| Cell metadata | Zenodo 3625024 | https://zenodo.org/record/3625024 |
| hs37d5 reference | HipSci FTP | https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/reference/hs37d5.fa.gz |
| GENCODE v19 GTF | GENCODE | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz |
| Imputed genotype BCFs | HipSci FTP | https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/data/vep_openaccess_bcf/ |
| WGS gVCF manifest | HipSci FTP | https://ftp.hipsci.ebi.ac.uk/hipsci/ftp/archive_datasets/ENA.ERP017015.wgs.normals.analysis_files.tsv |
| HipSci WGS gVCFs | ENA PRJEB15299 | https://www.ebi.ac.uk/ena/browser/view/PRJEB15299 |
| GitHub pipeline | singlecell_endodiff_paper | https://github.com/single-cell-genetics/singlecell_endodiff_paper |

---

## What to do about the 21 managed-access donors

To access cicb, dixh, eoxi, fasu, fejf, guss, guyj, hegp, koqx, nocf, nudd, oebj, oibg, oojs, pulk, tavh, tout, ueah, veku, walu, zoio you need:

1. Apply to the HipSci Data Access Committee via EGA: https://ega-archive.org/studies/EGAS00001001465
2. Once approved, download their scRNA CRAMs from EGA dataset EGAD00001005741
3. Their WGS VCFs are at EGA study EGAS00001001465 (same DAC)

The application asks for: your institution, intended use, IRB/ethics reference. Turnaround is typically 2–4 weeks.
