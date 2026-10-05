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


LOO_FILES = {"edgeR": "edger.tsv", "DESeq2": "deseq2.tsv", "limma": "limma.tsv"}
LOO_VERIFY_TOL = 1e-6


def _read_de(path: str) -> pd.DataFrame:
    t = pd.read_csv(path, sep="\t")
    if "log2FC" not in t.columns and "log2FoldChange" in t.columns:
        t = t.rename(columns={"log2FoldChange": "log2FC"})
    t["circ_id"] = t["circ_id"].astype(str)
    return t


def _sig_mask(t: pd.DataFrame, de_sig_by: str, fdr: float, lfc: float) -> pd.Series:
    col, thr = _resolve_sig(t, de_sig_by, fdr)
    return t[col].notna() & (t[col] < thr) & (t["log2FC"].abs() > lfc)


def _effect_stability(cands_ids: list, loo_dir: str, published: dict, min_ratio: float) -> pd.DataFrame:
    """Per-candidate, per-method effect-size stability under leave-one-sample-out.

    unstable (per method) = any dropped-sample run with sign(l) != sign(full), or |l| < min_ratio*|full|.
    Runs where the circRNA is absent (not testable) are counted in absent_<m> and do not set flags.
    Before any flag is used, the NONE (full-data) LOO run must reproduce the published log2FC of
    every candidate within LOO_VERIFY_TOL; otherwise this raises instead of ranking on a wrong baseline.
    """
    base = Path(loo_dir)
    if not (base / "NONE").is_dir():
        sys.exit(f"[rank] ERROR: --loo-dir {loo_dir} has no NONE run")
    drops = sorted(p.name for p in base.iterdir() if p.is_dir() and p.name != "NONE")
    if not drops:
        sys.exit(f"[rank] ERROR: --loo-dir {loo_dir} has no dropped-sample runs")
    out = pd.DataFrame(index=pd.Index(cands_ids, name="circ_id"))
    for m, fname in LOO_FILES.items():
        full = pd.read_csv(base / "NONE" / fname, sep="\t").set_index("circ_id")
        absent_full = [c for c in cands_ids if c not in full.index]
        if absent_full:
            sys.exit(f"[rank] ERROR: {len(absent_full)} candidates absent from {m} full LOO run "
                     f"(e.g. {absent_full[:3]}); refusing to rank")
        f = full.loc[cands_ids, "log2FC"].astype(float)
        pub_tab = published[m].set_index("circ_id")
        pub = pub_tab["log2FC"].reindex(cands_ids).astype(float)
        pub_p = pub_tab["pvalue"].reindex(cands_ids).astype(float)
        p_full = full.loc[cands_ids, "pvalue"].astype(float)
        dev = (f - pub).abs()
        dev_p = (p_full - pub_p).abs()
        # Tolerance gate on BOTH log2FC and p-value; the actual maxima are always printed.
        # Exact equality is not required (DESeq2 reproduces to ~1e-9).
        print(f"[rank] LOO gate {m}: max|dlog2FC|={dev.max():.3g} max|dp|={dev_p.max():.3g} "
              f"(tol {LOO_VERIFY_TOL:g})", file=sys.stderr)
        if dev.isna().any() or dev_p.isna().any() \
                or (dev > LOO_VERIFY_TOL).any() or (dev_p > LOO_VERIFY_TOL).any():
            n_bad = int(((dev > LOO_VERIFY_TOL) | (dev_p > LOO_VERIFY_TOL) | dev.isna() | dev_p.isna()).sum())
            sys.exit(f"[rank] ERROR: {m} full LOO run deviates from published results for {n_bad} candidates "
                     f"(max|dlog2FC|={dev.max():.3g}, max|dp|={dev_p.max():.3g}, tol {LOO_VERIFY_TOL:g}); "
                     f"refusing to rank")
        ratios, flips, shrinks, absent = [], pd.Series(0, index=cands_ids), pd.Series(0, index=cands_ids), pd.Series(0, index=cands_ids)
        for d in drops:
            t = pd.read_csv(base / d / fname, sep="\t").set_index("circ_id")
            l = t["log2FC"].reindex(cands_ids).astype(float)
            present = l.notna()
            absent += (~present).astype(int)
            flips += (present & (l * f < 0)).astype(int)
            r = l.abs() / f.abs()
            shrinks += (present & (r < min_ratio)).astype(int)
            ratios.append(r.rename(d))
        min_r = pd.concat(ratios, axis=1).min(axis=1, skipna=True)
        out[f"minratio_{m}"] = min_r.values
        out[f"sign_flip_{m}"] = flips.values
        out[f"drop_gt_{m}"] = shrinks.values
        out[f"absent_{m}"] = absent.values
        out[f"stable_{m}"] = ((flips == 0) & (shrinks == 0)).values
    out["effect_stable"] = out[[f"stable_{m}" for m in LOO_FILES]].all(axis=1)
    return out.reset_index()


def _write_run_params_filter(output: str, primary_method: str, cross_method_min: int, loo_dir: str,
                             stability_min_ratio: float, stability_max_candidates: int) -> None:
    """Record the biomarker filter settings in <results>/run_params.yaml (under biomarker_filter)."""
    run_params_path = Path(output).parent.parent / "run_params.yaml"
    existing = {}
    if run_params_path.exists():
        existing = yaml.safe_load(run_params_path.read_text(encoding="utf-8")) or {}
    existing["biomarker_filter"] = {
        "primary_method": primary_method,
        "cross_method_min": int(cross_method_min),
        "stability_min_ratio": float(stability_min_ratio),
        "stability_max_candidates": int(stability_max_candidates),
        "loo_dir": loo_dir or None,
        "output": Path(output).name,
    }
    run_params_path.write_text(yaml.safe_dump(existing, sort_keys=False, allow_unicode=True), encoding="utf-8")


def _write_status(path, status, n_candidates, n_three, n_pass, stability_run, loo_dir,
                  cross_method_min, stability_min_ratio) -> None:
    """Explicit run status (JSON). Written whenever --status-json is given, including the
    empty-result cases, so an empty table is never mistaken for a successful run."""
    if not path:
        return
    import json as _json
    doc = {
        "status": status,
        "n_candidates": int(n_candidates),
        "n_three_method_sig": int(n_three),
        "n_pass_default": int(n_pass),
        "stability_run": bool(stability_run),
        "loo_dir": loo_dir or None,
        "cross_method_min": int(cross_method_min),
        "stability_min_ratio": float(stability_min_ratio),
    }
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(_json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
    print(f"[rank] status={status} → {path}", file=sys.stderr)


def rank_biomarkers(
    de_file:           str,
    annot_file:        str,
    output:            str,
    fdr:               float = 0.05,
    lfc:               float = 1.0,
    de_sig_by:         str   = "pvalue",
    interactions_file: str   = "",
    de_deseq2:         str   = "",
    de_limma:          str   = "",
    sig_cap:           float = 10.0,
    fc_cap:            float = 5.0,
    component_weights: dict  = None,
    fully_evaluated_output: str = None,
    loo_dir:           str   = "",
    primary_method:    str   = "edgeR",
    de_edger:          str   = "",
    isoforms_file:     str   = "",
    stability_min_ratio: float = 0.5,
    cross_method_min:  int   = 3,
    write_run_params:  bool  = True,
    status_json:       str   = None,
    stability_max_candidates: int = 500,
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

    # Host gene / region from the isoform table (assign_isoforms.py output). Written into the
    # candidate file so the report does not have to join it itself.
    if isoforms_file:
        iso = pd.read_csv(isoforms_file, sep="\t", dtype=str,
                          usecols=["circ_id", "gene_name", "region"]).drop_duplicates("circ_id")
        merged = (merged.drop(columns=[c for c in ("gene_name", "region") if c in merged.columns])
                        .merge(iso, on="circ_id", how="left"))
        merged["gene_name"] = merged["gene_name"].fillna("")
        merged["region"] = merged["region"].fillna("")

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

    # No significant circRNA: write a HEADER-ONLY table plus an explicit status, and
    # exit 0 only as a report of that status (never write the full DE table as candidates).
    if sig.empty:
        print("[rank] WARNING: no significant circRNAs; writing header-only output.", file=sys.stderr)
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        merged.iloc[0:0].to_csv(output, sep="\t", index=False)
        _write_status(status_json, "no_significant_circrna", 0, 0, 0, bool(loo_dir), loo_dir,
                      cross_method_min, stability_min_ratio)
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
    # Strict: an unreadable or missing other-method table is an error, not a warning,
    # because the cross-method flag below must not silently fall back to fewer methods.
    # All three method tables are used symmetrically. `de` is the table of the PRIMARY
    # method (the one whose score is shown); the other two come from their own files.
    # The cross-method count and the LOO gate do not depend on which method is primary.
    path_of = {"edgeR": de_edger if primary_method != "edgeR" else de_file,
               "DESeq2": de_deseq2 if primary_method != "DESeq2" else de_file,
               "limma": de_limma if primary_method != "limma" else de_file}
    tables = {primary_method: de}
    for name in ("edgeR", "DESeq2", "limma"):
        if name == primary_method or not path_of[name]:
            continue
        if not Path(path_of[name]).exists():
            sys.exit(f"[rank] ERROR: {name} DE table given but missing: {path_of[name]}")
        tables[name] = _read_de(path_of[name])
    n_methods_available = len(tables)
    if cross_method_min > n_methods_available:
        sys.exit(f"[rank] ERROR: --cross-method-min {cross_method_min} but only "
                 f"{n_methods_available} DE method table(s) supplied (need --de-edger/--de-deseq2/--de-limma)")
    if loo_dir and n_methods_available < 3:
        sys.exit("[rank] ERROR: --loo-dir requires all three DE method tables (LOO gate checks each method)")
    sig_sets = {name: set(t.loc[_sig_mask(t, de_sig_by, fdr, lfc), "circ_id"].astype(str))
                for name, t in tables.items()}
    sig["n_sig_methods"] = sig["circ_id"].astype(str).apply(
        lambda c: sum(c in s for s in sig_sets.values())).astype(int)
    published = tables
    sig["three_method_sig"] = sig["n_sig_methods"] >= cross_method_min
    counts = sig["n_sig_methods"].value_counts().sort_index()
    for n, c in counts.items():
        print(f"[rank]   {c} circRNAs significant in {n}/3 method(s)", file=sys.stderr)

    # ── Effect-size stability (filter first, then test) ──────────────────────
    # Scores are normalised over the full significant set and are NOT recomputed on
    # the survivors. LOO is evaluated only for three-method candidates. The LOO fits
    # themselves are dataset-wide (one fit per dropped sample), so their cost does not
    # depend on how many candidates are tested.
    sig["effect_stable"] = pd.NA
    stability_run = bool(loo_dir)
    if stability_run:
        pub_tables = {k: v.assign(circ_id=v["circ_id"].astype(str)) for k, v in published.items()}
        survivors = sig.loc[sig["three_method_sig"], "circ_id"].astype(str).tolist()
        if len(survivors) > stability_max_candidates:
            sys.exit(f"[rank] ERROR: {len(survivors)} three-method candidates exceed "
                     f"--stability-max-candidates {stability_max_candidates}; raise the cap deliberately "
                     f"(stability is evaluated per candidate; the LOO fits themselves are per dataset)")
        if survivors:
            stab = _effect_stability(survivors, loo_dir, pub_tables, stability_min_ratio)
            sig["circ_id"] = sig["circ_id"].astype(str)
            sig = sig.drop(columns=["effect_stable"]).merge(stab, on="circ_id", how="left")
            sig["effect_stable"] = sig["effect_stable"].astype("boolean")
        else:
            sig["circ_id"] = sig["circ_id"].astype(str)
        print(f"[rank] stability: {len(survivors)} three-method candidates tested "
              f"(min ratio {stability_min_ratio}) from {loo_dir}", file=sys.stderr)
        sig["pass_default"] = sig["three_method_sig"] & sig["effect_stable"].fillna(False).astype(bool)
    else:
        sig["pass_default"] = sig["three_method_sig"]
        print("[rank] stability NOT run (no --loo-dir): pass_default = three-method only",
              file=sys.stderr)
    sig["pass_default"] = sig["pass_default"].astype(bool)

    # Deterministic order: score descending, ties broken by circ_id ascending.
    # (Plain sort_values is not stable on ties, so equal-score ranks could change between runs.)
    sig = sig.sort_values(["biomarker_score", "circ_id"], ascending=[False, True])
    sig.insert(sig.columns.get_loc("biomarker_score") + 1, "rank",
               range(1, len(sig) + 1))
    passing = sig["pass_default"]
    sig["rank_default"] = pd.NA
    sig.loc[passing, "rank_default"] = list(range(1, int(passing.sum()) + 1))

    Path(output).parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(output, sep="\t", index=False)
    n_pass = int(passing.sum())
    status = "ok" if n_pass > 0 else "no_candidates_pass"
    if n_pass == 0:
        print(f"[rank] WARNING: 0 of {len(sig)} candidates pass the default filter "
              f"(status=no_candidates_pass); all candidates are still written", file=sys.stderr)
    _write_status(status_json, status, len(sig), int(sig["three_method_sig"].sum()), n_pass,
                  stability_run, loo_dir, cross_method_min, stability_min_ratio)
    print(
        f"[rank] {len(sig)} biomarker candidates → {output}\n"
        f"[rank] top score: {sig['biomarker_score'].iloc[0]:.3f}  "
        f"({sig['circ_id'].iloc[0] if 'circ_id' in sig.columns else ''})",
        file=sys.stderr,
    )

    if write_run_params:
        _write_interaction_coverage(output, sig)
        _write_run_params_filter(output, primary_method, cross_method_min, loo_dir,
                                 stability_min_ratio, stability_max_candidates)

    if fully_evaluated_output:
        # Evidence-level-consistent ranking: only candidates where every
        # enabled component (including mirna/rbp) was actually evaluated,
        # not renormalized around a missing value. Scores are identical to
        # the main output (queried rows are unaffected by excluding
        # unqueried ones -- see CLAUDE.md diagnostic); only `rank` changes
        # because the comparison pool is smaller.
        fe = sig[sig["interaction_queried"]].copy()
        fe = fe.sort_values(["biomarker_score", "circ_id"], ascending=[False, True])
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
    parser.add_argument("--de-sig-by",  default="pvalue",
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
    parser.add_argument("--loo-dir", default="",
                        help="Directory from scripts/loo_stability.R (NONE/ + one dir per dropped "
                             "sample). If given, effect-size stability is tested for three-method "
                             "candidates; otherwise effect_stable is NA.")
    parser.add_argument("--stability-min-ratio", type=float, default=0.5,
                        help="A dropped-sample run fails if |log2FC_run| < ratio*|log2FC_full| "
                             "(default 0.5, i.e. a drop of more than 50%%). Sign flips always fail.")
    parser.add_argument("--cross-method-min", type=int, default=3, choices=[1, 2, 3],
                        help="Minimum number of DE methods (of those supplied) in which a candidate "
                             "must be significant to count as cross-method (default 3 = all supplied).")
    parser.add_argument("--primary-method", default="edgeR", choices=["edgeR", "DESeq2", "limma"],
                        help="Which method's table is --de (the scored table). The cross-method "
                             "flag and LOO gate are the same whichever method is primary.")
    parser.add_argument("--de-edger", default="",
                        help="edgeR DE table (required unless edgeR is the primary, in which case --de is used)")
    parser.add_argument("--isoforms", default="",
                        help="isoform_groups.tsv (assign_isoforms.py): adds gene_name and region to the candidate file")
    parser.add_argument("--no-run-params", action="store_true",
                        help="Do not write interaction coverage into <results>/run_params.yaml. "
                             "Required for v2 runs so published results directories are not modified.")
    parser.add_argument("--stability-max-candidates", type=int, default=500,
                        help="Upper bound on three-method candidates tested for stability (default 500); "
                             "exceeding it is an error, not a silent truncation.")
    parser.add_argument("--status-json", default=None,
                        help="Write an explicit run status JSON (status, counts, parameters) here.")
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
        loo_dir            = args.loo_dir,
        primary_method     = args.primary_method,
        de_edger           = args.de_edger,
        isoforms_file      = args.isoforms,
        stability_min_ratio = args.stability_min_ratio,
        cross_method_min   = args.cross_method_min,
        write_run_params   = not args.no_run_params,
        status_json        = args.status_json,
        stability_max_candidates = args.stability_max_candidates,
        fully_evaluated_output = (
            args.fully_evaluated_output
            or (str(Path(args.output).with_name(
                    Path(args.output).stem + "_fully_evaluated" + Path(args.output).suffix))
                if args.fully_evaluated_only else None)
        ),
    )


if __name__ == "__main__":
    main()
