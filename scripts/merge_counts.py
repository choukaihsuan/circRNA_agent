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
    """Extract circRNA IDs and BSJ read counts from one CIRIquant GTF."""
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

            id_match = re.search(r'circ_id\s+"([^"]+)"', attributes)
            circ_id = id_match.group(1) if id_match else f"{chrom}:{start}|{end}:{strand}"

            records.append({"circ_id": circ_id, "sample": sample_name, "BSJ": bsj})

    return pd.DataFrame(records, columns=["circ_id", "sample", "BSJ"])


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
) -> pd.DataFrame:
    """Concatenate all per-sample data, optionally filter by BED, pivot to wide format."""
    frames = []
    for gtf_path in gtf_files:
        sample = Path(gtf_path).parent.name
        df = parse_gtf(gtf_path, sample)
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    long_df = pd.concat(frames, ignore_index=True)

    if bed_files:
        confident = load_high_confidence(bed_files)
        before = len(long_df["circ_id"].unique())
        long_df = long_df[long_df["circ_id"].isin(confident)]
        after = len(long_df["circ_id"].unique())
        print(f"[filter] High-confidence filter: {before} → {after} circRNAs", file=sys.stderr)

    matrix = long_df.pivot_table(
        index="circ_id",
        columns="sample",
        values="BSJ",
        aggfunc="sum",
        fill_value=0,
    )
    matrix.columns.name = None
    matrix.index.name = "circ_id"
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="Build BSJ count matrix from CIRIquant GTFs")
    parser.add_argument("--gtfs", nargs="+", required=True, help="CIRIquant GTF files")
    parser.add_argument("--output", required=True, help="Output TSV matrix path")
    parser.add_argument("--filter-bed", nargs="*", dest="filter_bed",
                        help="High-confidence BED files for filtering (optional)")
    args = parser.parse_args()

    matrix = build_matrix(args.gtfs, args.filter_bed or [])
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t")
    print(f"[OK] Count matrix ({matrix.shape[0]} circRNAs × {matrix.shape[1]} samples) → {args.output}")


if __name__ == "__main__":
    main()
