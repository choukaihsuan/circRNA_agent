"""Traceable source for the six components and score of the Table 2 candidates (GSE133998, edgeR).

Components are recomputed from the v2 candidate file with biomarker_components.py, using the same
normalisation populations as production (the whole significant set; M and R over queried rows only).
Each row's equal-weight score is checked against the stored biomarker_score (4 decimals) before
any value is written. The output keeps full precision; the paper rounds it.

usage: python3 table2_components.py <GSE133998_results/de/biomarker_candidates_v2_edgeR.tsv> <out_tsv>
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from biomarker_components import COMPONENTS, components, weighted_score  # noqa: E402

EQUAL = {c: 1.0 for c in COMPONENTS}
GATE_TOL = 1e-4


def main() -> None:
    src, out = sys.argv[1], sys.argv[2]
    cand = pd.read_csv(src, sep="\t")
    comp = components(cand)
    score = weighted_score(comp, EQUAL)
    dev = (score - cand["biomarker_score"].astype(float)).abs()
    if dev.isna().any() or (dev > GATE_TOL).any():
        sys.exit(f"[table2] GATE FAILED: max |score - stored| = {dev.max():.3g}")
    print(f"[table2] gate OK: {len(cand)} rows, max |score - stored| = {dev.max():.2e}", file=sys.stderr)

    sel = cand["pass_default"].astype(str).isin(["True", "true"])
    t = pd.DataFrame({
        "rank_default": cand.loc[sel, "rank_default"].astype(int),
        "circ_id": cand.loc[sel, "circ_id"],
        "gene_name_assigned": cand.loc[sel, "gene_name"].fillna(""),
        "circbase_gene": cand.loc[sel, "circbase_gene"].fillna("novel"),
        "circbase_id": cand.loc[sel, "circbase_id"].fillna("novel"),
        "log2FC": cand.loc[sel, "log2FC"].astype(float),
        "pvalue": cand.loc[sel, "pvalue"].astype(float),
        "S": comp.loc[sel, "S"], "F": comp.loc[sel, "F"], "C": comp.loc[sel, "C"],
        "K": comp.loc[sel, "K"], "M": comp.loc[sel, "M"], "R": comp.loc[sel, "R"],
        "score_recomputed": score.loc[sel],
        "score_stored": cand.loc[sel, "biomarker_score"].astype(float),
        "rank_among_significant": cand.loc[sel, "rank"].astype(int),
    }).sort_values("rank_default")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(out, sep="\t", index=False, float_format="%.8g")
    print(t.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
