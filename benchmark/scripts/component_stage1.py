"""Stage 1 of the component validation (no external API calls).

1a  coverage of C, K, M, R on the SRP011042 RNase R reference (5,176 circRNAs; labels ER>1.5 / ER<0.5)
    and the positive rate in each defined subset.
1b  dependence among the annotation components on the 15 dual-tool candidate files (queried rows):
    P(in_circBase | n_mirna>0), P(n_mirna>0 | in_circBase=0), the same for n_rbp, and the
    point-biserial correlations K-M, K-R.
1c  AUROC of C versus log2(BSJ+1) against the ER labels, on the same rows.
    AUROC is reported without forcing a direction (values below 0.5 are reported as they are).

usage (server): python3 component_stage1.py <ground_truth.tsv> <our_method_summary.tsv> <circbase.txt> \
                                             <results_root> <out_tsv>
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "scripts"))  # pipeline scripts/
from annotate_circbase import _annotate_row, load_circbase  # noqa: E402

GSES = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
        "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289", "SRP156355"]
SLOP = 10


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney AUROC with average ranks for ties. No direction is imposed."""
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = pd.Series(np.concatenate([pos, neg])).rank(method="average").values
    r_pos = ranks[: len(pos)].sum()
    return float((r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def auroc_ci(auc: float, n1: int, n0: int) -> tuple[float, float]:
    """Hanley & McNeil (1982) standard error, normal 95% interval (clipped to [0, 1])."""
    if not (0 < auc < 1) or n1 < 2 or n0 < 2:
        return float("nan"), float("nan")
    q1 = auc / (2 - auc)
    q2 = 2 * auc * auc / (1 + auc)
    se = math.sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc * auc) + (n0 - 1) * (q2 - auc * auc)) / (n1 * n0))
    return max(0.0, auc - 1.96 * se), min(1.0, auc + 1.96 * se)


def fuzzy_find(df: pd.DataFrame, chrom: str, start: int, end: int, col: str):
    sub = df[df["chr"] == chrom]
    if sub.empty:
        return np.nan
    d = np.maximum((sub["start"] - start).abs(), (sub["end"] - end).abs())
    i = d.idxmin()
    return sub.loc[i, col] if d[i] <= SLOP else np.nan


def part_1a_1c(gt_path: str, summary_path: str, circbase_path: str, rows: list) -> None:
    gt = pd.read_csv(gt_path, sep="\t")
    summ = pd.read_csv(summary_path, sep="\t")
    summ = summ.rename(columns={"start": "start", "end": "end"})
    cb = load_circbase(circbase_path)
    n = len(gt)
    lab = gt["is_true"].values
    labeled = lab != -1

    # K: circBase match within slop, via annotate_circbase._annotate_row (same rule as production)
    k_vals = []
    for r in gt.itertuples(index=False):
        cid, _ = _annotate_row(r.chr, int(r.start), int(r.end), cb, SLOP)
        k_vals.append(0 if cid == "novel" else 1)
    gt["K"] = k_vals

    # C: consensus confidence_score for the same locus (SRR444655 circDEX consensus), if present
    gt["C"] = [fuzzy_find(summ, r.chr, int(r.start), int(r.end), "confidence_score") for r in gt.itertuples(index=False)]
    gt["C_defined"] = gt["C"].notna()
    gt["log2BSJ"] = np.log2(gt["bsj_total"].astype(float) + 1)

    def coverage(name: str, mask: np.ndarray) -> None:
        m_lab = mask & labeled
        n_pos = int(((lab == 1) & m_lab).sum())
        n_neg = int(((lab == 0) & m_lab).sum())
        rate = n_pos / (n_pos + n_neg) if (n_pos + n_neg) else float("nan")
        rows.append({"section": "1a", "dataset": "SRP011042_reference", "item": f"{name}_defined_n",
                     "value": int(mask.sum()), "n": n})
        rows.append({"section": "1a", "dataset": "SRP011042_reference", "item": f"{name}_defined_labeled_n",
                     "value": int(m_lab.sum()), "n": n})
        rows.append({"section": "1a", "dataset": "SRP011042_reference", "item": f"{name}_positive_rate",
                     "value": round(rate, 4), "n": int(n_pos + n_neg)})

    coverage("all", np.ones(n, bool))
    coverage("K_eq_1", gt["K"].values == 1)
    coverage("K_eq_0", gt["K"].values == 0)
    coverage("C", gt["C_defined"].values)
    coverage("C_undefined", ~gt["C_defined"].values)
    rows.append({"section": "1a", "dataset": "SRP011042_reference", "item": "M_defined_n", "value": "not computed (stage 2)", "n": n})
    rows.append({"section": "1a", "dataset": "SRP011042_reference", "item": "R_defined_n", "value": "not computed (stage 2)", "n": n})

    # 1c: C vs log2(BSJ+1) on the rows where C is defined and labeled
    m = labeled & gt["C_defined"].values
    y = lab[m]
    c = gt.loc[m, "C"].astype(float).values
    b = gt.loc[m, "log2BSJ"].values
    a_c, a_b = auroc(c, y), auroc(b, y)
    lo_c, hi_c = auroc_ci(a_c, int((y == 1).sum()), int((y == 0).sum()))
    lo_b, hi_b = auroc_ci(a_b, int((y == 1).sum()), int((y == 0).sum()))
    rho = float(pd.Series(c).corr(pd.Series(b), method="spearman"))
    rows += [
        {"section": "1c", "dataset": "SRP011042_reference(C_defined)", "item": "AUROC_C", "value": round(a_c, 4), "n": int(m.sum()), "ci_low": round(lo_c, 4), "ci_high": round(hi_c, 4)},
        {"section": "1c", "dataset": "SRP011042_reference(C_defined)", "item": "AUROC_log2BSJ", "value": round(a_b, 4), "n": int(m.sum()), "ci_low": round(lo_b, 4), "ci_high": round(hi_b, 4)},
        {"section": "1c", "dataset": "SRP011042_reference(C_defined)", "item": "spearman_C_vs_log2BSJ", "value": round(rho, 4), "n": int(m.sum())},
    ]
    # full labeled set for log2(BSJ+1) alone (C not defined for most rows)
    ml = labeled
    a_full = auroc(gt.loc[ml, "log2BSJ"].values, lab[ml])
    lo, hi = auroc_ci(a_full, int((lab[ml] == 1).sum()), int((lab[ml] == 0).sum()))
    rows.append({"section": "1c", "dataset": "SRP011042_reference(all_labeled)", "item": "AUROC_log2BSJ",
                 "value": round(a_full, 4), "n": int(ml.sum()), "ci_low": round(lo, 4), "ci_high": round(hi, 4)})
    print(f"[1c] C-defined labeled n={int(m.sum())}: AUROC_C={a_c:.4f} [{lo_c:.3f},{hi_c:.3f}]  "
          f"AUROC_log2BSJ={a_b:.4f} [{lo_b:.3f},{hi_b:.3f}]  spearman(C,log2BSJ)={rho:.3f}", file=sys.stderr)


def part_1b(results_root: str, rows: list) -> None:
    frames = []
    for g in GSES:
        p = Path(results_root) / f"{g}_results" / "de" / "biomarker_candidates_v2_edgeR.tsv"
        d = pd.read_csv(p, sep="\t")
        d = d[d["interaction_queried"].astype(str).isin(["True", "true"])].copy()
        d["dataset"] = g
        frames.append(d[["dataset", "circ_id", "in_circbase", "n_mirna", "n_rbp"]])
    pooled = pd.concat(frames, ignore_index=True)
    pooled["K"] = pooled["in_circbase"].astype(float)
    pooled["M"] = pooled["n_mirna"].astype(float)
    pooled["R"] = pooled["n_rbp"].astype(float)
    for name, d in [("pooled_15", pooled)] + [(g, pooled[pooled.dataset == g]) for g in GSES]:
        k1 = d[d.K == 1]
        k0 = d[d.K == 0]
        p_m_k1 = float((k1.M > 0).mean()) if len(k1) else float("nan")
        p_m_k0 = float((k0.M > 0).mean()) if len(k0) else float("nan")
        p_r_k1 = float((k1.R > 0).mean()) if len(k1) else float("nan")
        p_r_k0 = float((k0.R > 0).mean()) if len(k0) else float("nan")
        n_m0_pos = int(((d.K == 0) & (d.M > 0)).sum())
        n_r0_pos = int(((d.K == 0) & (d.R > 0)).sum())
        k_in_m = float((d[d.M > 0].K).mean()) if (d.M > 0).any() else float("nan")
        k_in_r = float((d[d.R > 0].K).mean()) if (d.R > 0).any() else float("nan")
        r_km = float(d["K"].corr(d["M"])) if len(d) > 2 else float("nan")
        r_kr = float(d["K"].corr(d["R"])) if len(d) > 2 else float("nan")
        r_km_bin = float(d["K"].corr((d["M"] > 0).astype(float))) if len(d) > 2 else float("nan")
        r_kr_bin = float(d["K"].corr((d["R"] > 0).astype(float))) if len(d) > 2 else float("nan")
        items = {
            "n_queried": len(d), "P(in_circBase|n_mirna>0)": round(k_in_m, 4),
            "P(in_circBase|n_rbp>0)": round(k_in_r, 4),
            "P(n_mirna>0|in_circBase=1)": round(p_m_k1, 4), "P(n_mirna>0|in_circBase=0)": round(p_m_k0, 4),
            "P(n_rbp>0|in_circBase=1)": round(p_r_k1, 4), "P(n_rbp>0|in_circBase=0)": round(p_r_k0, 4),
            "count(in_circBase=0 & n_mirna>0)": n_m0_pos, "count(in_circBase=0 & n_rbp>0)": n_r0_pos,
            "pointbiserial_K_M": round(r_km, 4), "pointbiserial_K_R": round(r_kr, 4),
            "pointbiserial_K_Mbin": round(r_km_bin, 4), "pointbiserial_K_Rbin": round(r_kr_bin, 4),
        }
        for k, v in items.items():
            rows.append({"section": "1b", "dataset": name, "item": k, "value": v, "n": len(d)})
    print(pooled.groupby("dataset").size().to_string(), file=sys.stderr)


def main() -> None:
    gt, summ, cb, root, out = sys.argv[1:6]
    rows: list = []
    part_1a_1c(gt, summ, cb, rows)
    part_1b(root, rows)
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"[stage1] written {len(df)} rows → {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
