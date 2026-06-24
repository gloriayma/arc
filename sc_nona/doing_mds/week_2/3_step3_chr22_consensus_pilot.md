# Step 3 (chr22 pilot) — donor consensus FASTA via streamed gVCF

Goal: build joxm's personalised chr22 reference FASTA without first downloading the 11.7 GB whole-genome gVCF. The plan's spec was pseudo-haploid with hs37d5-compatible coordinates: apply HOM-ALT calls only, leave hets at REF, preserve hs37d5 positions for easy sanity checking.

This doc records what works, what doesn't, and the three filter gotchas that wasted iterations before the consensus matched expectations.

## Final pipeline (works, ~30 s end-to-end)

```fish
cd data/hipsci/joxm/chr22

URL='https://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ448/ERZ448202/HPSI0114pf-joxm.wgs.gatk.haplotype_caller.20161201.genotypes.vcf.gz'

# Stream chr22 → HOM-ALT only → biallelic 1bp SNVs only (preserves hs37d5 coords)
bcftools view -r 22 \
    -i 'GT="1/1" && QUAL>=30 && FORMAT/GQ>=20' \
    --trim-alt-alleles \
    "$URL" 2>/dev/null \
  | bcftools view -m2 -M2 -i 'strlen(REF)==1 && strlen(ALT)==1' \
                  -Oz -o joxm.chr22.homalt_snps.vcf.gz
bcftools index -t joxm.chr22.homalt_snps.vcf.gz

# Extract chr22 from hs37d5 (50 MB vs 3 GB whole genome — faster for iteration)
samtools faidx ../../../reference/hs37d5.fa 22 > hs37d5.chr22.fa
samtools faidx hs37d5.chr22.fa

# Apply variants
bcftools consensus -f hs37d5.chr22.fa joxm.chr22.homalt_snps.vcf.gz \
  > joxm.chr22.consensus_snps.fa
samtools faidx joxm.chr22.consensus_snps.fa
```

Output:

- `joxm.chr22.homalt_snps.vcf.gz` — 16,418 strict 1bp HOM-ALT SNVs, 361 KB
- `joxm.chr22.consensus_snps.fa` — 51,304,566 bp (identical length to hs37d5 chr22), with 16,418 substitutions baked in
- Spot-checks at rs6518413 (16052239 A→G), rs5748728 (17387757 G→A), rs2070467 (24452885 A→G) all agree with hs37d5 → joxm consensus → VCF

The streaming part takes ~30 s and downloads only the bytes the .tbi says are inside `region 22:1-` plus a little overhead.

## Indexed remote access from `bcftools`

`bcftools` (via htslib) speaks HTTPS directly. With `-r 22` it queries the remote `.tbi` index for byte offsets and only fetches the matching ranges. Even though the gVCF is 11.7 GB, the chr22 stream pulls a small subset.

Two practical notes:

- If the local `.tbi` is in a different path than the gVCF, bcftools will *also* fetch the remote `.tbi` to the CWD (saw a duplicate file named `HPSI0114pf-joxm.wgs.gatk.haplotype_caller.20161201.genotypes.vcf.gz.tbi` appear next to the output VCF — md5-identical to the local one, removed). To force using a local index: `URL##idx##/abs/path/to/local.tbi`.
- For the whole-genome calls extract you'd be streaming all 11.7 GB anyway, so resume the curl download at that point.

## The three filter gotchas that wasted time

### 1. `-f PASS` drops everything in a single-sample GATK gVCF

The FILTER column in `HPSI0114pf-joxm.wgs.vcf.gz` is `.` for every record — GATK gVCF mode doesn't set FILTER at HaplotypeCaller stage. Filtering happens later via VQSR or hard filters, but the per-sample gVCFs HipSci releases haven't been run through that step.

```
=== FILTER value distribution on chr22 (first 10K records) ===
10000 .
```

`bcftools view -f PASS` filters to records where FILTER == "PASS" exactly, so it discarded all 10K records. **Use `QUAL` and `FORMAT/GQ` thresholds instead** (e.g. `QUAL>=30 && FORMAT/GQ>=20`), or run VQSR / `bcftools filter -e ...` on the gVCF first.

### 2. `<NON_REF>` sentinel makes every gVCF record look multi-allelic

GATK gVCF records carry `<NON_REF>` as a placeholder allele alongside any real ALT, so the ALT column looks like `T,<NON_REF>` even for a simple T substitution. `bcftools view -m2 -M2` (biallelic only) sees three alleles (REF + 2 ALT) and drops the record.

Fix:

```
bcftools view --trim-alt-alleles ...
```

This removes ALT alleles not seen in any genotype, including `<NON_REF>`. A cleaner alternative for a production pipeline is to run `gatk GenotypeGVCFs` first to convert the gVCF into a proper VCF — but for streaming a single chromosome, `--trim-alt-alleles` is fine.

### 3. `-v snps` is not strict; it leaks length-mismatched records that shift coordinates

After getting the first two fixes right, the consensus FASTA was 51,304,546 bp — **20 bp shorter than hs37d5 chr22**. Spot-checking the third variant (rs2070467) revealed the joxm consensus sequence at that position was completely unrelated to hs37d5: coords had drifted.

Cause: `bcftools view -v snps` keeps records bcftools classifies as "SNP-like" but with length-mismatched REF/ALT:

```
=== REF/ALT length distribution under -v snps ===
16418 1 -> 1   ← true SNVs
    3 1 -> 2   ← one bp + an insertion
    1 1 -> 5
    1 1 -> 8
    4 2 -> 1   ← two bp + a deletion
    1 2 -> 2
    1 29 -> 1  ← 29 bp -> 1 bp = 28 bp deletion!
    1 3 -> 1
```

Twelve length-mismatched records caused the 20 bp net shift. The fix is to drop the `-v snps` flag and filter explicitly by character length:

```
bcftools view -i 'strlen(REF)==1 && strlen(ALT)==1' ...
```

After this, 16,418 strict 1 bp SNVs applied, joxm consensus chr22 length matches hs37d5 chr22 length exactly, and rs2070467 reads correctly.

## What's left in `data/hipsci/joxm/chr22/` after cleanup

```
hs37d5.chr22.fa                 50 MB  chr22 extract of hs37d5
hs37d5.chr22.fa.fai             20 B
joxm.chr22.homalt_snps.vcf.gz   361 KB 16,418 strict 1bp HOM-ALT SNVs (PASS replaced with QUAL>=30, GQ>=20)
joxm.chr22.homalt_snps.vcf.gz.tbi 14 KB
joxm.chr22.consensus_snps.fa    50 MB  joxm chr22 with 16,418 substitutions, hs37d5 coords preserved
joxm.chr22.consensus_snps.fa.fai 20 B
consensus_snps.log              95 B   "Applied 16430 variants" trace
```

## What this pilot validated

- Streaming WGS gVCFs from ENA over HTTPS works; the `.tbi` is small enough to download separately and bcftools handles the rest.
- The filter recipe (`GT="1/1" && QUAL>=30 && FORMAT/GQ>=20` + `--trim-alt-alleles` + `strlen(REF)==1 && strlen(ALT)==1` + `-m2 -M2`) gives a coordinate-stable pseudo-haploid VCF that `bcftools consensus` can apply cleanly.
- The downstream pipeline (STAR index, CRAM alignment) can proceed against `joxm.chr22.consensus_snps.fa` immediately — no need to wait for the full-genome gVCF download.

## What the pilot deliberately deferred

- **Het sites**: skipped entirely under the pseudo-haploid model. Eventually want to handle them properly (either with IUPAC ambiguity codes via `bcftools consensus -I`, or with a full diploid pipeline as in the `2.2_genomes.md` `vcf2diploid`/`g2gtools` recommendations).
- **HOM-ALT indels**: dropped, since they shift hs37d5 coordinates. For the diploid version we'd want to keep them and track a chain file via `bcftools consensus -c CHAIN`.
- **VQSR / multi-sample joint calling**: the per-sample gVCF has no filter tags, so we're trusting QUAL/GQ thresholds. For production-grade calls, you'd run joint genotyping across the WGS cohort and apply VQSR.

These are all fine for a single-cell-line single-chromosome overfit test.
