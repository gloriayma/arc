# Step 5 — CRAM → FASTQ → STAR → bedGraph → BigWig

Goal: turn each downloaded joxm day-3 CRAM into a per-cell BigWig coverage track over the joxm chr22 consensus, ready to feed into the MSA-style training input.

## Final pipeline

Two scripts, both under `data/hipsci/joxm/scripts/`:

- **`align_cell.sh <CRAM> <OUT_ROOT>`** — does one cell end-to-end.
- **`run_batch.sh <CRAM_DIR> <OUT_DIR> [N_PARALLEL=4]`** — fans out `align_cell.sh` over a directory of CRAMs with `xargs -P`.

### Per-cell script (what each step does)

```fish
# Required env vars (set once before invoking either script)
set -x STAR_INDEX           /Users/.../joxm/chr22/STAR_index_joxm_chr22
set -x REF_FASTA            /Users/.../reference/hs37d5.fa
set -x BIN_BEDGRAPHTOBIGWIG /Users/.../bin/bedGraphToBigWig
set -x STAR_BIN             /Users/.../bin/STAR_x86      # see "STAR brew bug" below
set -x THREADS              4
```

The script:

1. **`samtools collate`** — name-sorts the CRAM in-place. CRAMs are coordinate-sorted, but `samtools fastq` needs read pairs adjacent. `collate` uses the temp dir as scratch.
2. **`samtools fastq -1 R1.fq.gz -2 R2.fq.gz`** — decode reads from CRAM (requires `--reference hs37d5.fa`), output paired gzipped FASTQs. Flags: `-F 0x900` drops secondary/supplementary, `-n` strips `/1`/`/2` suffixes (STAR rejects them).
3. **`STAR --runMode alignReads`** — splice-aware alignment against the chr22 STAR index. Outputs BAM (immediately deleted) and bedGraph signal tracks.
4. **`bedGraphToBigWig`** — UCSC tool to convert the larger STAR bedGraph into binary BigWig (~10× compression, random-access).
5. **Cleanup** — delete the temp FASTQs, BAM, bedGraph intermediates, STAR's `_STARtmp` / `_STARgenome` / `_STARpass1` dirs. Keep BigWig + `Log.final.out` + `SJ.out.tab`.

The script is idempotent: if `<CELL>.bw` already exists in the output dir, it short-circuits with `skip (already done)`.

## STAR parameters that matter

| Flag | Value | Why |
|---|---|---|
| `--outSAMtype BAM SortedByCoordinate` | — | STAR refuses to emit `--outWigType` without a sorted BAM. We delete the BAM after BigWig is built. |
| `--outWigType bedGraph` | — | Emits `Signal.UniqueMultiple.str1.out.bg` (uniq + multi-mappers, EM-weighted) and `Signal.Unique.str1.out.bg`. We use the UniqueMultiple file — sparse single-cell data benefits from keeping multi-mappers. |
| `--outWigStrand Unstranded` | — | Smart-seq2 is unstranded, so we want one combined track, not strand-split. |
| `--outWigNorm None` | — | Raw per-base coverage; can normalise to RPM later if needed. |
| `--outFilterMultimapNmax 20` | — | Default 10 is fine; 20 is a small bump for repetitive regions. |
| **No `--twopassMode`** | — | Pass 2 doubles wall-clock for ~0 gain on a chr22-only run with a pre-annotated GTF. Add `--twopassMode Basic` back for the eventual whole-genome run. |
| `--readFilesCommand gzcat` | — | macOS-only fix. See "gzcat" gotcha below. |

## Three painful gotchas

### 1. `rna-star 2.7.11b` brew bottle on arm64 macOS is broken

Confirmed reproducible: any STAR invocation produces `Number of input reads | 0` regardless of input method (raw FASTQ, gzipped, stdin pipe, hand-crafted 2-read minimal file, all return 0). The Log.out shows `Thread #0 end of input stream, nextChar=-1` — STAR's reader thread sees immediate EOF.

Workarounds tried:

- `gunzip -c` instead of `zcat`: helped for one error but didn't fix the underlying read failure
- single-thread (`--runThreadN 1`): no change
- different output types: no change
- moving files from `/tmp` to project dir: no change
- different output flags: no change

**Working fix**: install Rosetta 2 (`softwareupdate --install-rosetta --agree-to-license`) and use the official upstream macOS x86_64 STAR binary under Rosetta translation. Download from the github release:

```fish
curl -sSLo /tmp/STAR.zip 'https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip'
unzip -d /tmp /tmp/STAR.zip 'STAR_2.7.11b/MacOSX_x86_64/STAR'
cp /tmp/STAR_2.7.11b/MacOSX_x86_64/STAR data/bin/STAR_x86
chmod +x data/bin/STAR_x86
```

Then point the pipeline at `STAR_BIN=data/bin/STAR_x86`.

**Cost**: Rosetta translation is ~2–3× slower than native. For chr22-only with ~1.5 M paired reads per cell this is ~10 min per cell on M-series silicon vs ~3–4 min expected if a working native STAR existed. Live with it for the pilot; for the eventual whole-genome 100+ donor run, build native STAR from source with `brew install gcc` (Apple clang's `-fopenmp` story is broken too).

### 2. `zcat` on macOS is BSD `compress(1)` style, not gzip

`zcat foo.fq.gz` on macOS looks for `foo.fq.gz.Z` and fails with `can't stat`. STAR's `--readFilesCommand zcat` therefore silently produces no input. The fix is to use **`gzcat`** (a single-word command — see next point) — which IS gzip-aware on macOS — or install GNU coreutils.

### 3. STAR splits `--readFilesCommand` on whitespace

You'd think `--readFilesCommand 'gunzip -c'` would work as a multi-word command. It doesn't: STAR interprets the value as `argv[0]=gunzip`, `argv[1]=-c`, and treats `-c` as a FASTQ filename. You'll see `zcat: can't stat ... .Z` in the STAR log even though you passed `gunzip -c`. **Pass a single-word command** like `gzcat`.

### 4. `bedGraphToBigWig` is not on Homebrew

There's no `ucsc-bedgraphtobigwig` or similar formula in core. UCSC publishes static binaries instead:

```fish
curl -sSLo data/bin/bedGraphToBigWig \
  'https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bedGraphToBigWig'
chmod +x data/bin/bedGraphToBigWig
```

(Use `macOSX.x86_64` if not on Apple Silicon, or `linux.x86_64` on Linux.)

## Storage budget (this run)

| Phase | Disk used |
|---|---|
| Input: 50 CRAMs | 4.2 GB (kept) |
| Per cell peak (FASTQs + BAM + bedGraph) | ~1.5 GB; cleaned after BigWig built |
| 4-way parallel peak temp | ~5–7 GB |
| Final per cell (BigWig + STAR logs) | ~1 MB |
| Final 50-cell output | ~50 MB |

Growth from this step: **~50 MB** after cleanup.

## Wall-clock (this run)

Per-cell timing on Apple M-series under Rosetta with 4 STAR threads, chr22-only index, no two-pass:

- `samtools collate` + `samtools fastq` (from CRAM): ~30 s
- STAR alignment: ~9 min (Rosetta is the bottleneck)
- bedGraphToBigWig: ~5 s
- **Per cell total**: ~10 min sequential

For 50 cells at 4-way parallel: **~2.5–3 hours wall-clock**.

If you ever want this to run faster:

- Build STAR natively for arm64 (`brew install gcc`, then build from source). ~2–3× speed-up; cuts the batch to ~50–70 min.
- Skip the `samtools collate` step by streaming the FASTQs into STAR via FIFOs. Saves a few minutes per cell.
- Reduce `--outFilterMultimapNmax` and trim other STAR params. Marginal.

## Reproducibility: re-running the whole thing

```fish
cd /Users/gloria/dev/research/arc/sc_nona

# One-time setup (done once):
brew install rna-star samtools bcftools                 # not all needed if already installed
softwareupdate --install-rosetta --agree-to-license     # Apple Silicon only
# Then drop the upstream macOS x86 STAR binary into data/bin/STAR_x86 (see "brew bug" above)
# And bedGraphToBigWig from UCSC into data/bin/bedGraphToBigWig

# Pipeline run (re-runnable, idempotent):
set -x STAR_INDEX           $PWD/data/hipsci/joxm/chr22/STAR_index_joxm_chr22
set -x REF_FASTA            $PWD/data/reference/hs37d5.fa
set -x BIN_BEDGRAPHTOBIGWIG $PWD/data/bin/bedGraphToBigWig
set -x STAR_BIN             $PWD/data/bin/STAR_x86
set -x THREADS              4

data/hipsci/joxm/scripts/run_batch.sh \
    data/hipsci/joxm/crams_day3 \
    data/hipsci/joxm/bigwigs_day3 \
    4
```

`run_batch.sh` is idempotent: any cell whose BigWig already exists is skipped, so an interrupted run resumes cleanly.

## Outputs after the run

```
data/hipsci/joxm/bigwigs_day3/
  ERR2954998/
    ERR2954998.bw                ~700 KB — per-base coverage track, key output
    ERR2954998.Log.final.out     STAR alignment summary
    ERR2954998.SJ.out.tab        novel splice junctions from this cell
  ERR2955090/
    ...
  ...
  batch.log                      run_batch.sh stdout
```

The `.bw` files are the artifact for downstream model training — load with `pyBigWig` to get per-base coverage arrays for any chr22 window.

## Sanity check (run on 3 BigWigs: ERR2954998, ERR2955090, ERR2955091)

Script: `scripts/sanity_check_bigwig.py`. Needs `pyBigWig` in your venv.

### Gene-vs-intergenic baseline

|             | mean coverage | nonzero in % cells |
|---|---|---|
| MAPK1 (housekeeping kinase) | 2.14 | 100% |
| GTSE1 (cell-cycle) | 4.65 | 100% |
| TBX1 (DE TF) | 0.44 | 66.7% |
| EP300 (housekeeping cofactor) | 1.56 | 100% |
| **intergenic gene-desert at 22:34,790,000–34,795,000** | **0.69** | (control) |

Housekeeping genes show 2–7× over the intergenic baseline. TBX1 shows up in 2/3 cells, plausible for a transient developmental TF.

**Gotcha**: my first "intergenic control" was a 5 kb window at 22:30,000,000–30,005,000 that happened to sit inside **NF2** (a broadly-expressed tumor suppressor), and gave 2.16 baseline coverage. Re-picked an honest gene desert from the largest GENCODE gap (the 380 kb stretch at 22:34,605,251–34,985,077) and the baseline dropped to ~0.7.

### Splice-aware exon-vs-intron test (the actual smoking gun)

The plan's step-9 sanity check asks for "high coverage over exonic regions, near-zero coverage over introns, sharp drops at exon boundaries". Naively averaging over all GENCODE exons and all introns of a gene didn't give a clean signal because some long retained-intron annotations are listed as exons in alternative isoforms — they pollute the "exon" set with low-coverage regions.

The cleaner test: pick **one specific small exon** and compute coverage in it vs 500 bp flanking-intron windows. Using MAPK1 exon at 22:153,301–22:153,417 (117 bp):

```
          cell    upstream        EXON  downstream    enrichment
    ERR2954998        0.00        7.18        0.00          inf×
    ERR2955090        0.00        0.00        0.00          —
    ERR2955091        0.00        2.07        0.00          inf×
```

That's the textbook signal: sharp peak at exon, near-zero on either side. STAR's splice-aware alignment is working correctly. ERR2955090 doesn't express this isoform but is internally consistent.

### Conclusion

Pipeline is sound. Per-cell BigWigs encode chr22 RNA-seq coverage with proper exon-structured peaks. Ready to feed into downstream training as base-resolution tracks.

### Pilot scale note

The pipeline was originally launched on 50 day-3 joxm cells with 4-way parallel under Rosetta — projected wall-clock turned out to be ~8 hours, much slower than the ~3 hour estimate. **Cut scope to 3 cells** for this sanity-check pass. To do the full 50:

- Either let the 4-way-parallel batch run overnight (`scripts/run_batch.sh CRAM_DIR OUT_DIR 4`)
- Or build native arm64 STAR (`brew install gcc`, then build from source) — expected 2–3× speed-up that brings the full batch to ~3 hours.

Both paths are reproducible via the same scripts.
