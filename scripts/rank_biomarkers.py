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
import yaml


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


DEFAULT_COMPONENT_WEIGHTS = {
    "sig":   {"enabled": True, "weight": 1.0},
    "fc":    {"enabled": True, "weight": 1.0},
    "conf":  {"enabled": True, "weight": 1.0},
    "known": {"enabled": True, "weight": 1.0},
    "mirna": {"enabled": True, "weight": 1.0},
    "rbp":   {"enabled": True, "weight": 1.0},
}
# known_bonus is a raw 0/1 flag, never min-max normalized (unchanged from
# the original 4D/6D formula) -- every other component is.
_UNNORMALIZED_COMPONENTS = {"known"}


def _compute_weighted_score(sig: pd.DataFrame, sig_col: str, sig_cap: float,
                             fc_cap: float, component_weights: dict = None) -> pd.Series:
    """
    Two-layer weighted biomarker score.

    Layer 1 (config-level): weights are normalized within the set of
    *enabled* components (`weight / sum(enabled weights)`), so equal weights
    on all six components reproduces a simple average regardless of how many
    are enabled -- this is this paper's default setting.

    Layer 2 (per-row): mirna/rbp are only defined for candidates where
    `interaction_queried` is True (see caller). For a row where they are
    missing, the row's weights are renormalized over whichever enabled
    components it actually has -- a missing value is excluded from the
    average, never silently treated as 0.
    """
    components = component_weights or DEFAULT_COMPONENT_WEIGHTS
    enabled = {name: float(cfg.get("weight", 1.0)) for name, cfg in components.items()
               if cfg.get("enabled", True) and float(cfg.get("weight", 1.0)) > 0}
    if not enabled:
        sys.exit("[rank] ERROR: no biomarker_score components are enabled with weight > 0 "
                  "-- check de.biomarker.components in config.yaml")

    raw = pd.DataFrame(index=sig.index)
    raw["sig"]   = sig[sig_col].apply(lambda p: min(-math.log10(max(p, 1e-10)), sig_cap))
    raw["fc"]    = sig["log2FC"].abs().clip(upper=fc_cap)
    raw["conf"]  = sig["confidence_score"].clip(lower=0)
    raw["known"] = sig["in_circbase"].astype(float)
    raw["mirna"] = sig["n_mirna"].where(sig["interaction_queried"])  # NaN where not queried
    raw["rbp"]   = sig["n_rbp"].where(sig["interaction_queried"])

    norm = pd.DataFrame(index=sig.index, columns=raw.columns, dtype=float)
    for c in raw.columns:
        if c in _UNNORMALIZED_COMPONENTS:
            norm[c] = raw[c]
            continue
        valid = raw[c].notna()
        if valid.any():
            norm.loc[valid, c] = _minmax(raw.loc[valid, c])
        # rows with no valid value in this component stay NaN -> excluded below

    weight_row = pd.Series({c: enabled.get(c, 0.0) for c in norm.columns})
    effective_weight = norm.notna().mul(weight_row, axis=1)
    row_weight_sum = effective_weight.sum(axis=1)
    if (row_weight_sum <= 0).any():
        n_bad = int((row_weight_sum <= 0).sum())
        print(f"[rank] WARNING: {n_bad} candidate(s) have no available enabled component "
              f"(all NaN) -- their biomarker_score will be NaN", file=sys.stderr)
        row_weight_sum = row_weight_sum.replace(0, pd.NA)

    weighted_sum = norm.fillna(0).mul(effective_weight).sum(axis=1)
    score = (weighted_sum / row_weight_sum).round(4)

    dims_str = " + ".join(sorted(enabled))
    n_unqueried = int((~sig["interaction_queried"]).sum())
    extra = ""
    if n_unqueried and ("mirna" in enabled or "rbp" in enabled):
        extra = (f"; {n_unqueried}/{len(sig)} candidates were never in the interaction "
                 f"query union (predict_interactions top-N) -- their score is renormalized "
                 f"over the remaining enabled components ('M/R not evaluated')")
    print(f"[rank] score dimensions (enabled, weight-normalized): {dims_str}{extra}",
          file=sys.stderr)
    return score


def _write_interaction_coverage(output_path: str, sig: pd.DataFrame) -> None:
    """
    Record miRNA/RBP query coverage into {results_dir}/run_params.yaml so
    report.html's Analysis parameters section (and anyone auditing the run
    later) can see what fraction of candidates the M/R components were
    actually evaluated for, without re-deriving it from biomarker_candidates.tsv.
    Merges into the existing file (if any) rather than overwriting it, since
    this is one of several sections that file is expected to accumulate.
    """
    if "interaction_queried" not in sig.columns or len(sig) == 0:
        return
    try:
        results_dir = Path(output_path).resolve().parent.parent
        run_params_path = results_dir / "run_params.yaml"
        existing: dict = {}
        if run_params_path.exists():
            try:
                existing = yaml.safe_load(run_params_path.read_text()) or {}
            except Exception as e:
                print(f"[rank] WARNING: could not parse existing {run_params_path} "
                      f"({e}); it will be overwritten", file=sys.stderr)
                existing = {}
        n_sig = int(len(sig))
        n_queried = int(sig["interaction_queried"].sum())
        existing["biomarker_interaction_coverage"] = {
            "n_significant": n_sig,
            "n_queried": n_queried,
            "coverage_pct": round(n_queried / n_sig * 100, 1) if n_sig else 0.0,
            "note": ("miRNA/RBP components of the biomarker_score are only evaluated for "
                     "candidates predict_interactions.py actually queried (top-N union); "
                     "candidates outside that union have their score renormalized over the "
                     "remaining enabled components, not penalized to 0."),
        }
        run_params_path.write_text(yaml.safe_dump(existing, sort_keys=False, allow_unicode=True))
        print(f"[rank] interaction coverage {n_queried}/{n_sig} "
              f"({existing['biomarker_interaction_coverage']['coverage_pct']}%) "
              f"-> {run_params_path}", file=sys.stderr)
    except OSError as e:
        print(f"[rank] WARNING: could not write run_params.yaml ({e}); "
              f"coverage info only available in biomarker_candidates.tsv", file=sys.stderr)


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
    sig_cap:           float = 10.0,
    fc_cap:            float = 5.0,
    component_weights: dict  = None,
    fully_evaluated_output: str = None,
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
    # A circ_id has a key in interactions.json iff predict_interactions.py
    # actually queried it (it was in the top-N union AND metadata resolved);
    # that is true even when the query came back with zero binders. Absence
    # of the key means "never queried" -- distinct from "queried, found none".
    # This distinction is what the two-layer weighting below depends on; see
    # CLAUDE.md "六維生物標記分數 — M/R 缺值陷阱" for why conflating the two
    # (the pre-2026-10 behavior) silently re-ranks by query coverage instead
    # of biology once mirna/rbp weights are raised.
    icounts = _load_interaction_counts(interactions_file) if interactions_file else {}
    queried_ids = set(icounts.keys())
    if "circ_id" in sig.columns:
        sig["interaction_queried"] = sig["circ_id"].isin(queried_ids)
        sig["n_mirna"] = sig["circ_id"].map(lambda x: icounts.get(x, (0, 0))[0])
        sig["n_rbp"]   = sig["circ_id"].map(lambda x: icounts.get(x, (0, 0))[1])
    else:
        sig["interaction_queried"] = False
        sig["n_mirna"] = 0
        sig["n_rbp"]   = 0

    sig["biomarker_score"] = _compute_weighted_score(sig, sig_col, sig_cap, fc_cap, component_weights)

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

    _write_interaction_coverage(output, sig)

    if fully_evaluated_output:
        # Evidence-level-consistent ranking: only candidates where every
        # enabled component (including mirna/rbp) was actually evaluated,
        # not renormalized around a missing value. Scores are identical to
        # the main output (queried rows are unaffected by excluding
        # unqueried ones -- see CLAUDE.md diagnostic); only `rank` changes
        # because the comparison pool is smaller.
        fe = sig[sig["interaction_queried"]].copy()
        fe = fe.sort_values("biomarker_score", ascending=False)
        fe["rank"] = range(1, len(fe) + 1)
        Path(fully_evaluated_output).parent.mkdir(parents=True, exist_ok=True)
        fe.to_csv(fully_evaluated_output, sep="\t", index=False)
        print(
            f"[rank] {len(fe)}/{len(sig)} fully-evaluated candidates (interaction_queried=True) "
            f"→ {fully_evaluated_output}",
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
    parser.add_argument("--sig-cap", type=float, default=10.0,
                        help="Cap for -log10(p) in sig_norm component (default: 10)")
    parser.add_argument("--fc-cap",  type=float, default=5.0,
                        help="Cap for |log2FC| in fc_norm component (default: 5)")
    parser.add_argument("--fully-evaluated-only", action="store_true",
                        help="Also write <output>_fully_evaluated.tsv (or --fully-evaluated-output "
                             "path) containing only candidates with interaction_queried=True, "
                             "re-ranked within that evidence-consistent subset. Does not replace "
                             "the main --output file.")
    parser.add_argument("--fully-evaluated-output", default=None,
                        help="Output path for the fully-evaluated-only TSV "
                             "(default: <output> with _fully_evaluated inserted before the extension)")
    parser.add_argument("--components-json", default=None,
                        help="JSON dict: {component: {enabled, weight}} for sig/fc/conf/"
                             "known/mirna/rbp. Default: all enabled, weight 1.0 (equal, "
                             "this paper's setting).")
    args = parser.parse_args()

    component_weights = None
    if args.components_json:
        try:
            component_weights = json.loads(args.components_json)
        except json.JSONDecodeError as e:
            sys.exit(f"[rank] ERROR: --components-json is not valid JSON: {e}")

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
        sig_cap            = args.sig_cap,
        fc_cap             = args.fc_cap,
        component_weights  = component_weights,
        fully_evaluated_output = (
            args.fully_evaluated_output
            or (str(Path(args.output).with_name(
                    Path(args.output).stem + "_fully_evaluated" + Path(args.output).suffix))
                if args.fully_evaluated_only else None)
        ),
    )


if __name__ == "__main__":
    main()
