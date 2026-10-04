"""Weight sensitivity of the six-component biomarker score (edgeR primary, 15 dual-tool datasets).

No new data and no new comparison: everything is recomputed from the existing v2 candidate files
(biomarker_candidates_v2_edgeR.tsv) with biomarker_components.py.

Gate first: the 'equal' scheme must reproduce the production biomarker_score column of the same
file row by row (within the 4-decimal rounding that rank_biomarkers writes). If any dataset fails
the gate, the script stops and writes nothing for the remaining schemes.

Schemes:
  equal        all six components weight 1 (production default)
  drop_S..drop_R  one component removed (weight 0), the others weight 1
  de_heavy     S,F weight 2; C,K,M,R weight 1
  anno_heavy   S,F weight 1; C,K,M,R weight 2
  de_only      S,F only (C,K,M,R weight 0)

Score rule: for each row, the weighted mean over the components that are valid (non-NaN) for that
row, with the denominator equal to the sum of the weights of those components. The denominator is
not a fixed sum of six. M and R are NaN for rows that were never interaction-queried.

Metrics per dataset and scheme (reference = equal):
  spearman_rho    Spearman correlation of scores across all candidates
  top10_overlap   how many of the equal top-10 are in the scheme's top-10
  top1_changed    whether the top-ranked candidate differs
  n_candidates    denominator for interpretation

usage: python3 weight_sensitivity.py <results_root> <out_tsv> <summary_tsv>
  results_root holds <GSE>_results/de/biomarker_candidates_v2_edgeR.tsv (one per dataset)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from biomarker_components import COMPONENTS, components, weighted_score  # noqa: E402

GSES = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
        "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289", "SRP156355"]
GATE_TOL = 1e-4  # production writes the score rounded to 4 decimals

EQUAL = {c: 1.0 for c in COMPONENTS}
SCHEMES = {"equal": EQUAL}
for c in COMPONENTS:
    SCHEMES[f"drop_{c}"] = {k: (0.0 if k == c else 1.0) for k in COMPONENTS}
SCHEMES["de_heavy"] = {"S": 2, "F": 2, "C": 1, "K": 1, "M": 1, "R": 1}
SCHEMES["anno_heavy"] = {"S": 1, "F": 1, "C": 2, "K": 2, "M": 2, "R": 2}
SCHEMES["de_only"] = {"S": 1, "F": 1, "C": 0, "K": 0, "M": 0, "R": 0}


def top_ids(score: pd.Series, ids: pd.Series, k: int = 10) -> list:
    order = pd.DataFrame({"s": score, "id": ids}).sort_values(["s", "id"], ascending=[False, True])
    return order["id"].head(k).tolist()


def main() -> None:
    root, out_tsv, out_sum = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = []
    for gse in GSES:
        path = Path(root) / f"{gse}_results" / "de" / "biomarker_candidates_v2_edgeR.tsv"
        cand = pd.read_csv(path, sep="\t")
        comp = components(cand)
        eq = weighted_score(comp, EQUAL)
        prod = cand["biomarker_score"].astype(float)
        dev = (eq - prod).abs()
        if dev.isna().any() or (dev > GATE_TOL).any():
            sys.exit(f"[weights] GATE FAILED for {gse}: max |equal - production| = {dev.max():.3g} "
                     f"(tol {GATE_TOL:g}); nothing further written")
        print(f"[weights] gate {gse}: n={len(cand)} max|equal-production|={dev.max():.2e}", file=sys.stderr)
        eq_top = top_ids(eq, cand["circ_id"])
        for name, w in SCHEMES.items():
            sc = weighted_score(comp, w)
            rho = eq.corr(sc, method="spearman")
            sc_top = top_ids(sc, cand["circ_id"])
            rows.append({
                "dataset": gse, "method": "edgeR_ciriquant", "scheme": name,
                "n_candidates": int(len(cand)),
                "spearman_rho": round(float(rho), 4),
                "top10_overlap": len(set(eq_top) & set(sc_top)),
                "top1_changed": bool(eq_top[0] != sc_top[0]),
                "M_missing_fraction": round(float(comp["M"].isna().mean()), 4),
            })
    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep="\t", index=False)

    summ = []
    for name in SCHEMES:
        sub = df[df.scheme == name]
        summ.append({"scheme": name,
                     "rho_median": round(sub.spearman_rho.median(), 4),
                     "rho_min": round(sub.spearman_rho.min(), 4),
                     "rho_max": round(sub.spearman_rho.max(), 4),
                     "top10_median": sub.top10_overlap.median(),
                     "top10_min": int(sub.top10_overlap.min()),
                     "top10_max": int(sub.top10_overlap.max()),
                     "top1_changed_n": int(sub.top1_changed.sum()),
                     "datasets_rho_below_0.9": ";".join(sub.loc[sub.spearman_rho < 0.9, "dataset"])})
    pd.DataFrame(summ).to_csv(out_sum, sep="\t", index=False)
    print(pd.DataFrame(summ).to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
