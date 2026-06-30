"""Emit one BED6 file per cell from a scRNA gene-count CSV.

Score column = expression count for that gene in that cell.
Genes with zero counts are skipped.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path


def load_gene_coords(gtf_path: Path) -> dict[str, tuple[str, int, int, str]]:
    """Map ENSG id (version-stripped) -> (chrom, start_1based, end_1based, strand)."""
    coords: dict[str, tuple[str, int, int, str]] = {}
    opener = gzip.open if gtf_path.suffix == ".gz" else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, strand, attrs = f[0], int(f[3]), int(f[4]), f[6], f[8]
            gid = None
            for a in attrs.split(";"):
                a = a.strip()
                if a.startswith("gene_id "):
                    gid = a.split('"')[1].split(".")[0]
                    break
            if gid is not None:
                coords[gid] = (chrom, start, end, strand)
    return coords


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True, type=Path)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--num-cells", type=int, default=3,
                    help="Number of cells (columns) to export. Use -1 for all.")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading gene coordinates from {args.gtf} ...", file=sys.stderr)
    coords = load_gene_coords(args.gtf)
    print(f"  {len(coords)} gene records", file=sys.stderr)

    with open(args.counts, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        cell_ids = header[1:]
        if args.num_cells > 0:
            cell_ids = cell_ids[: args.num_cells]
        n_cells = len(cell_ids)
        print(f"Exporting {n_cells} cell(s): {cell_ids}", file=sys.stderr)

        records: list[list[tuple[str, int, int, str, float, str]]] = [[] for _ in range(n_cells)]
        n_missing = 0
        n_genes = 0
        for row in reader:
            n_genes += 1
            gene_field = row[0]
            ensg = gene_field.split("_", 1)[0]
            name = gene_field.split("_", 1)[1] if "_" in gene_field else ensg
            hit = coords.get(ensg)
            if hit is None:
                n_missing += 1
                continue
            chrom, start, end, strand = hit
            for i in range(n_cells):
                v = float(row[1 + i])
                if v == 0:
                    continue
                # BED is 0-based half-open
                records[i].append((chrom, start - 1, end, name, v, strand))

    print(f"  {n_genes} genes read, {n_missing} unmatched to GTF", file=sys.stderr)

    for cid, recs in zip(cell_ids, records):
        safe = cid.replace("#", "_").replace("/", "_")
        out = args.outdir / f"cell_{safe}.bed"
        recs.sort(key=lambda r: (r[0], r[1]))
        with open(out, "w") as fh:
            fh.write(f'track name="{cid}" description="scRNA counts for cell {cid}" '
                     f'useScore=0 visibility=full\n')
            for chrom, s, e, name, v, strand in recs:
                fh.write(f"{chrom}\t{s}\t{e}\t{name}\t{v:g}\t{strand}\n")
        print(f"  wrote {out} ({len(recs)} non-zero genes)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
