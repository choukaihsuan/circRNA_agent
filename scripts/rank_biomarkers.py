"""
rank_biomarkers.py – Rank DE circRNAs by composite biomarker score.

Composite score combines up to six evidence dimensions:

    1. significance  = -log10(p-value),  capped at 10
    2. effect_size   = |log2FC|,          capped at 5
    3. confidence    = confidence_score from consensus_filter (0–∞, normalised 0–1)
    4. known_bonus   = 1 if in circBase, 0 otherwise
    5. n_mirna       = distinct predicted miRNA binders (from interactions.json, if provided)
    6. n_rbp         = distinct predicted RBP binders  (from interactions.json, if provided)

    biomarker_score = (sig_norm + fc_norm + conf_norm + known_bonus [+ mirna_n + rbp_n]) / N
    where N = 4 without interactions, 6 with interactions.

Each normalised dimension is min-max scaled within the significant set.
The known_bonus is not normalised (0 or 1). CircRNAs absent from interactions.json
receive n_mirna = 0 and n_rbp = 0.

Output TSV columns (all input DE columns + added):
    confidence_score, in_circbase, circbase_id, circbase_gene,
    n_mirna, n_rbp, biomarker_score, rank

Usage:
    python scripts/rank_biomarkers.py \
        --de           de/de_results.tsv \
        --annot        circRNA/circbase_annotated.tsv \
        --output       de/biomarker_candidates.tsv \
        --interactions de/interactions.json \
        --fdr     0.05 \
        --lfc     1.0
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import pandas as pd


def _minmax(series: pd.Series) -> pd.Series:
    lo, hi = series.min(), series.max()
    if hi == lo:
        return pd.Series(0.5, index=series.index)
    return (series - lo) / (hi - lo)


def _is_in_circ(entry: dict, key: str) -> bool:
    """True if this interaction entry actually binds within the circRNA exon sequence."""
    if entry.get("source") == "CircInteractome":
        return True   # CI entries are always circRNA-specific
    return bool(entry.get("in_circ", False))


def _load_interaction_counts(interactions_file: str) -> dict:
    """Return {circ_id: (n_mirna, n_rbp)} counting only in-circRNA binding sites."""
    try:
        with open(interactions_file) as f:
            data = json.load(f)
    except Exception:
        return {}
    counts = {}
    for circ_id, val in data.items():
        mirnas = val.get("mirna") or []
        rbps   = val.get("rbp")   or []
        n_mirna = len(set(
            m.get("miRNAName", "") for m in mirnas
            if m.get("miRNAName") and _is_in_circ(m, "miRNAName")
        ))
        n_rbp = len(set(
            r.get("RBPName", "") for r in rbps
            if r.get("RBPName") and _is_in_circ(r, "RBPName")
        ))
        counts[circ_id] = (n_mirna, n_rbp)
    return counts


def _resolve_sig(de: pd.DataFrame, de_sig_by: str, fdr: float):
    """Return (column, threshold) based on cascade logic matching analysis.R."""
    if de_sig_by == "auto":
        if "qvalue" in de.columns and (de["qvalue"] < 0.2).any():
            return "qvalue", 0.2
        return "pvalue", 0.05
    if de_sig_by == "qvalue":
        return "qvalue", 0.2
    if de_sig_by == "pvalue":
        return "pvalue", fdr
    return "padj", fdr


def rank_biomarkers(
    de_file:           str,
    annot_file:        str,
    output:            str,
    fdr:               float = 0.05,
    lfc:               float = 1.0,
    de_sig_by:         str   = "auto",
    interactions_file: str   = "",
    de_deseq2:         str   = "",
    de_limma:          str   = "",
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
    sig_col, sig_thr = _resolve_sig(merged, de_sig_by, fdr)
    if sig_col not in merged.columns or "log2FC" not in merged.columns:
        sys.exit(f"[rank] ERROR: DE file must contain {sig_col} and log2FC columns")
    print(f"[rank] significance: {sig_col} < {sig_thr}", file=sys.stderr)

    sig_mask = (
        merged[sig_col].notna() &
        (merged[sig_col] < sig_thr) &
        (merged["log2FC"].abs() > lfc)
    )
    sig = merged.loc[sig_mask].copy()

    if sig.empty:
        print("[rank] No significant circRNAs; writing empty output.", file=sys.stderr)
        merged.to_csv(output, sep="\t", index=False)
        return

    # ── Interaction counts (optional) ─────────────────────────────────────────
    icounts = _load_interaction_counts(interactions_file) if interactions_file else {}
    if icounts and "circ_id" in sig.columns:
        sig["n_mirna"] = sig["circ_id"].map(lambda x: icounts.get(x, (0, 0))[0])
        sig["n_rbp"]   = sig["circ_id"].map(lambda x: icounts.get(x, (0, 0))[1])
    else:
        sig["n_mirna"] = 0
        sig["n_rbp"]   = 0

    use_interactions = icounts and sig["n_mirna"].max() + sig["n_rbp"].max() > 0

    # ── Component scores ───────────────────────────────────────────────────────
    sig["_sig"] = sig[sig_col].apply(lambda p: min(-math.log10(max(p, 1e-10)), 10))
    sig["_fc"]  = sig["log2FC"].abs().clip(upper=5)
    sig["_conf"] = sig["confidence_score"].clip(lower=0)

    sig["_sig_n"]  = _minmax(sig["_sig"])
    sig["_fc_n"]   = _minmax(sig["_fc"])
    sig["_conf_n"] = _minmax(sig["_conf"])

    if use_interactions:
        sig["_mirna_n"] = _minmax(sig["n_mirna"])
        sig["_rbp_n"]   = _minmax(sig["n_rbp"])
        sig["biomarker_score"] = (
            sig["_sig_n"] + sig["_fc_n"] + sig["_conf_n"] + sig["in_circbase"]
            + sig["_mirna_n"] + sig["_rbp_n"]
        ) / 6.0
        sig = sig.drop(columns=["_sig", "_fc", "_conf", "_sig_n", "_fc_n", "_conf_n",
                                 "_mirna_n", "_rbp_n"])
        print("[rank] score dimensions: sig + fc + conf + circbase + mirna + rbp (6D)",
              file=sys.stderr)
    else:
        sig["biomarker_score"] = (
            sig["_sig_n"] + sig["_fc_n"] + sig["_conf_n"] + sig["in_circbase"]
        ) / 4.0
        sig = sig.drop(columns=["_sig", "_fc", "_conf", "_sig_n", "_fc_n", "_conf_n"])
        print("[rank] score dimensions: sig + fc + conf + circbase (4D)", file=sys.stderr)

    sig["biomarker_score"] = sig["biomarker_score"].round(4)

    # ── n_sig_methods: count across all three DE methods ──────────────────────
    sig["n_sig_methods"] = 1  # primary method is always 1 (these are already sig)
    for other_file in [de_deseq2, de_limma]:
        if not other_file or not Path(other_file).exists():
            continue
        try:
            other = pd.read_csv(other_file, sep="\t")
            if "log2FC" not in other.columns and "log2FoldChange" in other.columns:
                other = other.rename(columns={"log2FoldChange": "log2FC"})
            o_col, o_thr = _resolve_sig(other, de_sig_by, fdr)
            if o_col in other.columns and "log2FC" in other.columns and "circ_id" in other.columns:
                other_ids = set(other.loc[
                    other[o_col].notna() &
                    (other[o_col] < o_thr) &
                    (other["log2FC"].abs() > lfc), "circ_id"
                ].tolist())
                sig.loc[sig["circ_id"].isin(other_ids), "n_sig_methods"] += 1
        except Exception as e:
            print(f"[rank] warning: could not read {other_file}: {e}", file=sys.stderr)

    sig["n_sig_methods"] = sig["n_sig_methods"].astype(int)
    counts = sig["n_sig_methods"].value_counts().sort_index()
    for n, c in counts.items():
        print(f"[rank]   {c} circRNAs significant in {n}/3 method(s)", file=sys.stderr)

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
    parser.add_argument("--de",           required=True, help="DE results TSV (from analysis.R)")
    parser.add_argument("--annot",        required=True, help="circBase-annotated summary TSV")
    parser.add_argument("--output",       required=True, help="Output ranked TSV")
    parser.add_argument("--interactions", default="",
                        help="interactions.json from predict_interactions (optional; adds miRNA/RBP dimensions)")
    parser.add_argument("--de-deseq2",   default="", help="DESeq2 DE results TSV (optional; for n_sig_methods)")
    parser.add_argument("--de-limma",    default="", help="limma DE results TSV (optional; for n_sig_methods)")
    parser.add_argument("--fdr",          type=float, default=0.05)
    parser.add_argument("--lfc",          type=float, default=1.0)
    parser.add_argument("--de-sig-by",  default="auto",
                        choices=["auto", "pvalue", "padj", "qvalue"],
                        help="Significance column: auto=cascade (q<0.2→p<0.05), "
                             "pvalue/padj/qvalue=fixed")
    # backward compat: --use-pvalue maps to --de-sig-by pvalue
    parser.add_argument("--use-pvalue", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()

    de_sig_by = args.de_sig_by
    if args.use_pvalue and de_sig_by == "auto":
        de_sig_by = "pvalue"

    rank_biomarkers(
        de_file            = args.de,
        annot_file         = args.annot,
        output             = args.output,
        fdr                = args.fdr,
        lfc                = args.lfc,
        de_sig_by          = de_sig_by,
        interactions_file  = args.interactions,
        de_deseq2          = args.de_deseq2,
        de_limma           = args.de_limma,
    )


if __name__ == "__main__":
    main()
