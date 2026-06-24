# Step 1 — metadata download for the single-donor pilot

Completed step 1 of `own_writing/8_bigger_picture_plan.md`. Picked donor **`joxm`** (HipSci line `HPSI0114i-joxm_1`) as the pilot:

- Highest cell count among open-access Cuomo donors: **1,415 cells** in Zenodo, **1,328 with matching ENA runs** (the other 87 are in the EGA managed-access half).
- Has a parent-fibroblast WGS gVCF (`HPSI0114pf-joxm`, ERZ448202), which is the germline genome — cleaner than the iPSC clone for building a personalised reference.
- Cells span all four days (day0/1/2/3) so will be useful for the eventual covariation experiments.

---

## Files produced

```
data/cuomo/
  cell_metadata_cols.tsv          # 39 MB — Zenodo dump, all 36,044 post-QC cells × 90+ columns
  cell_to_donor.tsv               # cell_id → donor_stem / donor_long_id / day, all 36,044 cells
  ena_runs.tsv                    # 46,468 runs from PRJEB14362 with CRAM + FASTQ URLs
  cell_to_ena.tsv                 # 23,793 Zenodo cells joined to their ENA run + CRAM/FASTQ URLs
  cell_to_ena_joxm.tsv            # 1,328 joxm cells with CRAM + FASTQ URLs
  cell_to_ena_missing.tsv         # 12,251 Zenodo cells with NO ENA match (EGA managed-access half)
  build_donor_cram_manifest.py    # builder script

data/hipsci/
  wgs_analysis_manifest.tsv       # 530 WGS per-sample gVCF analyses from PRJEB15299
  open_access_donors.tsv          # 125 Cuomo donors with status + preferred VCF URL
  build_open_access_table.py      # builder script
  joxm/
    HPSI0114pf-joxm.wgs.vcf.gz.tbi   # 2.7 MB tabix index (complete)
    HPSI0114pf-joxm.wgs.vcf.gz       # 3.4 GB / 11.7 GB partial — download paused
    joxm_vcf_download.log
```

The full-genome gVCF download was paused after 3.4 GB. Resume any time with:

```fish
cd data/hipsci/joxm
curl -L -C - -o HPSI0114pf-joxm.wgs.vcf.gz \
  'https://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ448/ERZ448202/HPSI0114pf-joxm.wgs.gatk.haplotype_caller.20161201.genotypes.vcf.gz'
```

For the chr22-only pilot you can also stream straight from ENA — see step 3 doc.

Three small builders sit alongside the data, so reproducing or extending is straightforward.

---

## Open-access summary (verified against `2.2b_access_verification.md`)

| Cuomo donor status | Count |
|---|---|
| Open WGS gVCF (parent fibroblast `pf`) | 101 |
| Open WGS gVCF (iPSC clone `i`, donor `giju` only) | 1 |
| Open imputed BCF only (`iezw`, `uenn`) | 2 |
| **Open-access total** | **104** |
| EGA managed-access (DAC application required) | 21 |

---

## Cell → ENA mapping rule (validated)

`cell_id "<run>_<lane>#<tag>"  ==  experiment_alias[len("SC_EXP_"):]`

The `#` is a Sanger tag character that must be URL-encoded as `%23` when fetching the CRAM (the manifest already has this encoded).

For joxm, 1,328 / 1,415 Zenodo cells have matching open ENA runs. The 87 missing cells all share Sanger run prefix `21843`, which is presumably in the EGA managed half — globally across all 36,044 Zenodo cells, **12,251 (34%)** are similarly EGA-only.

---

## joxm cell distribution

| Day | n cells (ENA-matched) |
|---|---|
| day0 | 383 |
| day1 | 453 |
| day2 | 317 |
| day3 | 175 |
| **Total** | **1,328** |

The plan calls for ~50 cells for the first overfit run; joxm gives plenty of choice across timepoints.

---

## Example CRAM URL (verified HTTP 200)

```
https://ftp.sra.ebi.ac.uk/vol1/run/ERR295/ERR2953552/24722_3%23101.cram
```

Content-Length: 27,344,784 bytes (~27 MB per cell). For all 1,328 joxm cells this is roughly **~40 GB of CRAMs** — comfortably feasible.

---

## What's next (step 2)

1. Download the `hs37d5.fa.gz` reference + GENCODE v19 GTF.
2. Decide whether to wait for the full WGS gVCF or stream chromosome 22 first for a quicker pipeline sanity check.

re's the command to stream chr22 from the remote VCF without downloading the full 11.7 GB:

 cd /Users/gloria/dev/research/arc/sc_nona/data/hipsci/joxm

      bcftools view -r 22 \
        https://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ448/ERZ448202/HPSI0114pf-joxm.wgs.gatk.haplotype_caller.20161201.genotypes.vcf.gz \
        -Oz -o joxm.chr22.vcf.gz

      bcftools index -t joxm.chr22.vcf.gz

    This uses the .tbi index (already downloaded at data/hipsci/joxm/HPSI0114pf-joxm.wgs.vcf.gz.tbi) to issue HTTP range requests and pull only the chr22 byte ranges from the remote file — should finish in a few minutes instead of ~45.

    Two preconditions worth knowing:
    1. bcftools must be installed. The access-verification doc noted it wasn't installed in the project venv. On macOS: brew install bcftools. Or via conda: conda install -c bioconda bcftools.
    2. The contig name is 22, not chr22 — hs37d5 uses bare numbers. If you pass chr22 you'll get an empty result.

    Optional follow-up to make it usable downstream:
    # Pre-filter to PASS biallelic non-HOM-REF calls (drops most of the gVCF bulk)
    bcftools view -f PASS -m2 -M2 -e 'GT="./." || GT="0/0"' \
      joxm.chr22.vcf.gz -Oz -o joxm.chr22.calls.vcf.gz
    bcftools index -t joxm.chr22.calls.vcf.gz

    After that you have everything you need to run bcftools consensus against the chr22 slice of hs37d5 and build the step-3 personalised FASTA for chr22 — without waiting for the full WGS download.


    
3. Begin step 3 (build the donor's reference FASTA via `bcftools consensus`).
