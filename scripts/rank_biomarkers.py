"""
rank_biomarkers.py – Rank DE circRNAs by composite biomarker score.

Composite score combines four evidence dimensions:

    1. significance  = -log10(padj),  capped at 10
    2. effect_size   = |log2FC|,       capped at 5
    3. confidence    = confidence_score from consensus_filter (0–∞, normalised 0–1)
    4. known_bonus   = 1 if in circBase, 0 otherwise

    biomarker_score = (sig_norm + fc_norm + conf_norm + known_bonus) / 4

Each dimension is min-max normalised within the significant set before summing,
so no single dimension dominates. The known_bonus is not normalised (0 or 1).

Output TSV columns (all input DE columns + added):
    confidence_score, in_circbase, circbase_id, circbase_gene, biomarker_score, rank

Usage:
    python scripts/rank_biomarkers.py \
        --de      de/de_results.tsv \
        --annot   circRNA/circbase_annotated.tsv \
        --output  de/biomarker_candidates.tsv \
        --fdr     0.05 \
        --lfc     1.0
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd


def _minmax(series: pd.Series) -> pd.Series:
    lo, hi = series.min(), series.max()
    if hi == lo:
        return pd.Series(0.5, index=series.index)
    return (series - lo) / (hi - lo)


def rank_biomarkers(
    de_file:    str,
    annot_file: str,
    output:     str,
    fdr:        float = 0.05,
    lfc:        float = 1.0,
) -> None:
    de    = pd.read_csv(de_file,    sep="\t")
    annot = pd.read_csv(annot_file, sep="\t")

    # Normalise log2FC column name (analysis.R always outputs log2FC)
    if "log2FC" not in de.columns and "log2FoldChange" in de.columns:
        de = de.rename(columns={"log2FoldChange": "log2FC"})

    # Merge annotation by circ_id (or chr/start/end fallback)
    annot_cols = [c for c in ["circ_id", "confidence_score", "in_circbase",
                               "circbase_id", "circbase_gene"] if c in annot.columns]
    if "circ_id" in annot.columns and "circ_id" in de.columns:
        merged = de.merge(annot[annot_cols], on="circ_id", how="left")
    else:
        merged = de.copy()
        for c in annot_cols:
            if c not in merged.columns:
                merged[c] = pd.NA

    # Fill missing annotation
    merged["in_circbase"]   = merged.get("in_circbase",   pd.Series(0, index=merged.index)).fillna(0).astype(int)
    merged["confidence_score"] = pd.to_numeric(merged.get("confidence_score", 0), errors="coerce").fillna(0)
    merged["circbase_id"]   = merged.get("circbase_id",   pd.Series("", index=merged.index)).fillna("")
    merged["circbase_gene"] = merged.get("circbase_gene", pd.Series("", index=merged.index)).fillna("")

    # Filter to significant DE circRNAs
    if "padj" not in merged.columns or "log2FC" not in merged.columns:
        sys.exit("[rank] ERROR: DE file must contain padj and log2FC columns")

    sig_mask = (
        merged["padj"].notna() &
        (merged["padj"] < fdr) &
        (merged["log2FC"].abs() > lfc)
    )
    sig = merged.loc[sig_mask].copy()

    if sig.empty:
        print("[rank] No significant circRNAs; writing empty output.", file=sys.stderr)
        merged.to_csv(output, sep="\t", index=False)
        return

    # ── Component scores ───────────────────────────────────────────────────────
    sig["_sig"] = sig["padj"].apply(lambda p: min(-math.log10(max(p, 1e-10)), 10))
    sig["_fc"]  = sig["log2FC"].abs().clip(upper=5)
    sig["_conf"] = sig["confidence_score"].clip(lower=0)

    sig["_sig_n"]  = _minmax(sig["_sig"])
    sig["_fc_n"]   = _minmax(sig["_fc"])
    sig["_conf_n"] = _minmax(sig["_conf"])

    sig["biomarker_score"] = (
        sig["_sig_n"] + sig["_fc_n"] + sig["_conf_n"] + sig["in_circbase"]
    ) / 4.0
    sig["biomarker_score"] = sig["biomarker_score"].round(4)

    sig = sig.drop(columns=["_sig", "_fc", "_conf", "_sig_n", "_fc_n", "_conf_n"])
    sig = sig.sort_values("biomarker_score", ascending=False)
    sig.insert(sig.columns.get_loc("biomarker_score") + 1, "rank",
               range(1, len(sig) + 1))

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(output, sep="\t", index=False)
    print(
        f"[rank] {len(sig)} biomarker candidates → {output}\n"
        f"[rank] top score: {sig['biomarker_score'].iloc[0]:.3f}  "
        f"({sig['circ_id'].iloc[0] if 'circ_id' in sig.columns else ''})",
        file=sys.stderr,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rank DE circRNAs by composite biomarker score"
    )
    parser.add_argument("--de",     required=True, help="DE results TSV (from analysis.R)")
    parser.add_argument("--annot",  required=True, help="circBase-annotated summary TSV")
    parser.add_argument("--output", required=True, help="Output ranked TSV")
    parser.add_argument("--fdr",    type=float, default=0.05)
    parser.add_argument("--lfc",    type=float, default=1.0)
    args = parser.parse_args()

    rank_biomarkers(
        de_file    = args.de,
        annot_file = args.annot,
        output     = args.output,
        fdr        = args.fdr,
        lfc        = args.lfc,
    )


if __name__ == "__main__":
    main()
