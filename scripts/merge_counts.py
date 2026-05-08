"""
merge_counts.py – Parse CIRIquant per-sample GTF outputs and build a
BSJ count matrix (circRNAs × samples).

Called as a Snakemake script: receives snakemake.input.gtfs and
snakemake.output.matrix via the snakemake object.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd


def parse_gtf(gtf_path: str, sample_name: str) -> pd.DataFrame:
    """
    Extract circRNA IDs and BSJ read counts from one CIRIquant GTF.

    CIRIquant GTF attribute fields include:
      circ_id "<id>"; ... BSJ <float>; ...
    """
    records = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "circRNA":
                continue

            chrom, _, _, start, end, _, strand, _, attributes = parts

            # Parse BSJ count from attributes
            bsj_match = re.search(r'BSJ\s+([\d.]+)', attributes)
            if not bsj_match:
                continue
            bsj = float(bsj_match.group(1))

            # Parse circ_id if present, else construct from coordinates
            id_match = re.search(r'circ_id\s+"([^"]+)"', attributes)
            circ_id = id_match.group(1) if id_match else f"{chrom}:{start}|{end}:{strand}"

            records.append({"circ_id": circ_id, "sample": sample_name, "BSJ": bsj})

    return pd.DataFrame(records, columns=["circ_id", "sample", "BSJ"])


def build_matrix(gtf_files: list[str]) -> pd.DataFrame:
    """Concatenate all per-sample data and pivot to wide format."""
    frames = []
    for gtf_path in gtf_files:
        sample = Path(gtf_path).parent.name
        df = parse_gtf(gtf_path, sample)
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    long_df = pd.concat(frames, ignore_index=True)
    matrix  = long_df.pivot_table(
        index="circ_id",
        columns="sample",
        values="BSJ",
        aggfunc="sum",
        fill_value=0,
    )
    matrix.columns.name = None
    matrix.index.name   = "circ_id"
    return matrix


# ── Snakemake entry point ────────────────────────────────────────────────────

gtf_files   = snakemake.input.gtfs   # type: ignore[name-defined]
output_file = snakemake.output.matrix  # type: ignore[name-defined]

matrix = build_matrix(list(gtf_files))
Path(output_file).parent.mkdir(parents=True, exist_ok=True)
matrix.to_csv(output_file, sep="\t")
print(f"[OK] Count matrix ({matrix.shape[0]} circRNAs × {matrix.shape[1]} samples) → {output_file}")
