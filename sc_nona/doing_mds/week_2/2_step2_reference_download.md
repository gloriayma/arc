# Step 2 — reference genome + GENCODE annotation

Goal: have the same `hs37d5` reference HipSci used to align their CRAMs, plus a GENCODE v19 annotation in the same coordinate system, so the chr22 sanity-check pipeline (step 3) can build a personalised consensus FASTA and STAR can splice-align cell reads against it later (step 3.5+).

## Files produced (all in `data/reference/`)

| File | Size | Source | Notes |
|---|---|---|---|
| `hs37d5.fa.gz` | 860 MB | `ftp.hipsci.ebi.ac.uk/hipsci/ftp/reference/hs37d5.fa.gz` | Kept as archive |
| `hs37d5.fa` | 3.0 GB | `gunzip -k` of the above | Working copy |
| `hs37d5.fa.fai` | 2.7 KB | `samtools faidx` | 86 contigs |
| `gencode.v19.annotation.gtf.gz` | 36 MB | `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19` | Original, `chr1`/`chrM` naming |
| `gencode.v19.annotation.nochr.gtf.gz` | 39 MB | derived | Remapped contig names — **use this with hs37d5** |
| `hs37d5_download.log` | 13 KB | curl progress | |

## Why hs37d5 specifically

HipSci aligned all their WGS CRAMs against this exact build (GRCh37 primary + decoys + EBV + rCRS mitochondrion). Using anything else (plain GRCh37, GRCh38, etc.) would either misalign reads or require liftover later. The reference path is in every GATK HaplotypeCaller header line of the per-sample gVCFs (`reference_sequence=/lustre/scratch116/vr/projects/hipsci/vrpipe/refs/human/ncbi37/hs37d5.fa`), confirming this is what they used.

## Contig naming — the silent compatibility trap

`hs37d5` uses **Ensembl-style** contig names: `1, 2, ..., 22, X, Y, MT, GL...`. GENCODE v19 ships with **UCSC-style** names: `chr1, chr2, ..., chr22, chrX, chrY, chrM`. The CRAMs from HipSci were encoded against hs37d5, so they use the no-prefix names too.

If you ever feed STAR a `chr*`-named GTF and a `1*`-named FASTA, alignment will produce zero spliced reads at known genes — STAR silently treats the chromosomes as different contigs. The remap rule is:

```
chr1..chr22 -> 1..22
chrX/chrY   -> X/Y
chrM        -> MT     ← the easy-to-miss swap
```

Done with one awk — **note the explicit `-F'\t'`**:

```fish
zcat < gencode.v19.annotation.gtf.gz | \
  awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} {
    sub(/^chr/, "", $1);
    if ($1 == "M") $1 = "MT";
    print
  }' | gzip > gencode.v19.annotation.nochr.gtf.gz
```

Verified the remapped GTF has exactly `{1..22, X, Y, MT}` for the data contigs.

**Don't omit the `-F'\t'`.** The GTF spec only uses tabs between the 9 mandatory fields; the 9th attribute field contains spaces (`gene_id "ENSG..."; transcript_id "..."`). Without `-F'\t'`, awk default-splits on any whitespace, modifies any field, then re-joins with `OFS="\t"` — silently replacing every space in the attribute string with a tab. Result: corrupted GTF, STAR will fail to parse it. This bit me on first build (see `3.5_step3.5_star_index.md`).

## What's in `hs37d5.fa`

86 contigs total in `.fai`:

- **25 standard**: 1–22, X, Y, MT
- **~60 decoy / unplaced contigs** (`GL000xxx.1`): assorted unplaced supercontigs and HuRef-derived decoys
- **EBV**: `NC_007605` (Epstein-Barr virus, type 1)
- **`hs37d5`**: the synthetic "decoy" itself (~35 Mb), built from BAC clones + NA12878

The decoys catch reads that would otherwise misalign to similar-but-non-canonical sequence in the genome; useful for variant calling. GENCODE has no annotations on the decoys (since they're not "the genome"), so the no-chr GTF only covers the 25 standard contigs.

## The "trailing garbage ignored" warning

```
gunzip: hs37d5.fa.gz: trailing garbage ignored
```

This is benign. `hs37d5.fa.gz` is a multi-stream gzip with some trailing padding bytes; gunzip prints the warning but extracts cleanly. Verified by checking the downloaded size matches the HEAD Content-Length exactly (901,788,727 bytes), and the .fai index lists all 86 contigs at the documented lengths (chr1 = 249,250,621 etc.).

## Why GENCODE v19 and not something newer

v19 is the last GENCODE release built against GRCh37. Anything from v20 onward is GRCh38-only and would not align coordinates to hs37d5. The Cuomo paper itself used v19. Practical implication for the project: when we later want to overlay annotations on the BigWig coverage tracks, the GTF coordinates already match the donor consensus FASTA (since the consensus preserves hs37d5 coordinates — see step-3 doc).

## Disk

- Reference + indexes: ~3.9 GB on disk
- Can drop `hs37d5.fa.gz` (860 MB) if disk pressure ever bites; everything downstream uses the uncompressed `.fa` + `.fai`.
