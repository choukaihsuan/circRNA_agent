"""
merge_counts.py – Parse CIRIquant per-sample GTF outputs and build a
BSJ count matrix (circRNAs × samples), optionally filtered by
high-confidence BED files from multi-tool consensus voting.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def parse_gtf(gtf_path: str, sample_name: str) -> pd.DataFrame:
    """Extract circRNA IDs, BSJ and FSJ read counts from one CIRIquant GTF."""
    records = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "circRNA":
                continue

            chrom, _, _, start, end, _, strand, _, attributes = parts

            bsj_match = re.search(r'BSJ\s+([\d.]+)', attributes)
            if not bsj_match:
                continue
            bsj = float(bsj_match.group(1))

            fsj_match = re.search(r'FSJ\s+([\d.]+)', attributes)
            fsj = float(fsj_match.group(1)) if fsj_match else 0.0

            id_match = re.search(r'circ_id\s+"([^"]+)"', attributes)
            circ_id = id_match.group(1) if id_match else f"{chrom}:{start}|{end}:{strand}"

            records.append({"circ_id": circ_id, "sample": sample_name, "BSJ": bsj, "FSJ": fsj})

    return pd.DataFrame(records, columns=["circ_id", "sample", "BSJ", "FSJ"])


def load_high_confidence(bed_files: list[str], slop: int = 0) -> set[str]:
    """Load circRNA IDs from BED files (chr:start-end format)."""
    ids: set[str] = set()
    for bed in bed_files:
        with open(bed) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.strip().split("\t")
                if len(cols) < 3:
                    continue
                circ_id = f"{cols[0]}:{cols[1]}|{cols[2]}"
                ids.add(circ_id)
    return ids


def build_matrix(
    gtf_files: list[str],
    bed_files: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (BSJ matrix, FSJ matrix), both circRNAs × samples."""
    frames = []
    for gtf_path in gtf_files:
        sample = Path(gtf_path).parent.name
        df = parse_gtf(gtf_path, sample)
        frames.append(df)

    if not frames:
        empty = pd.DataFrame()
        return empty, empty

    long_df = pd.concat(frames, ignore_index=True)

    if bed_files:
        confident = load_high_confidence(bed_files)
        before = len(long_df["circ_id"].unique())
        long_df = long_df[long_df["circ_id"].isin(confident)]
        after = len(long_df["circ_id"].unique())
        print(f"[filter] High-confidence filter: {before} → {after} circRNAs", file=sys.stderr)

    def _pivot(col: str) -> pd.DataFrame:
        m = long_df.pivot_table(
            index="circ_id", columns="sample", values=col,
            aggfunc="sum", fill_value=0,
        )
        m.columns.name = None
        m.index.name = "circ_id"
        return m

    return _pivot("BSJ"), _pivot("FSJ")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build BSJ/FSJ count matrices from CIRIquant GTFs")
    parser.add_argument("--gtfs", nargs="+", required=True, help="CIRIquant GTF files")
    parser.add_argument("--output",     required=True, help="Output BSJ count matrix TSV")
    parser.add_argument("--output-fsj", dest="output_fsj", default=None,
                        help="Output FSJ count matrix TSV (optional)")
    parser.add_argument("--filter-bed", nargs="*", dest="filter_bed",
                        help="High-confidence BED files for filtering (optional)")
    args = parser.parse_args()

    bsj_matrix, fsj_matrix = build_matrix(args.gtfs, args.filter_bed or [])
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    bsj_matrix.to_csv(args.output, sep="\t")
    print(f"[OK] BSJ matrix ({bsj_matrix.shape[0]} circRNAs × {bsj_matrix.shape[1]} samples) → {args.output}")

    if args.output_fsj:
        fsj_matrix.to_csv(args.output_fsj, sep="\t")
        print(f"[OK] FSJ matrix ({fsj_matrix.shape[0]} circRNAs × {fsj_matrix.shape[1]} samples) → {args.output_fsj}")


if __name__ == "__main__":
    main()
