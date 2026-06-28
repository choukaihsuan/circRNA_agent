"""
accuracy_benchmark.py – 偵測準確率比較（對照 RNase R ground truth）

Three multi-tool detection strategies evaluated (single-tool methods excluded):
  1. Our adaptive consensus filter
     CIRIquant + DCC, slop=10 bp, BSJ/FSJ pseudo-circ QC
  2. CirComPara2 simulation (Gaffo et al. 2022)
     CIRIquant + DCC, slop=10 bp, NO pseudo-circ QC
  3. nf-core/circrna simulation (Digby-Bell et al. 2023)
     CIRIquant + CIRCexplorer2 + find_circ, slop=0 (exact coords), min_tools=2

Metrics: Precision, Recall, F1, Specificity, AUC-PR
Stratification: BSJ count quartile in total RNA (Q1 / Q2-Q3 / Q4)
  Boundaries computed from ground-truth set (CirComPara2 / nf-core style)

Outputs:
  --output-summary    results/benchmark/accuracy_summary.tsv
  --output-stratified results/benchmark/stratified_f1.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── Coordinate matching ───────────────────────────────────────────────────────

def _parse_id(circ_id: str) -> tuple[str, int, int] | None:
    m = re.match(r'^(.+):(\d+)\|(\d+)$', circ_id)
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def _coord_match_in_set(circ_id: str, id_set: set[str], slop: int) -> bool:
    """Check if circ_id has a match within slop bp anywhere in id_set."""
    parsed = _parse_id(circ_id)
    if parsed is None:
        return circ_id in id_set
    chrom, start, end = parsed
    for k in id_set:
        kp = _parse_id(k)
        if kp is None or kp[0] != chrom:
            continue
        if max(abs(kp[1] - start), abs(kp[2] - end)) <= slop:
            return True
    return False


def _score_for_id(circ_id: str, scores: dict[str, float], id_set: set[str], slop: int) -> float:
    """Look up the score for circ_id via fuzzy matching."""
    parsed = _parse_id(circ_id)
    if parsed is None:
        return scores.get(circ_id, 0.0)
    chrom, start, end = parsed
    best_score, best_dist = 0.0, float("inf")
    for k, s in scores.items():
        kp = _parse_id(k)
        if kp is None or kp[0] != chrom:
            continue
        dist = max(abs(kp[1] - start), abs(kp[2] - end))
        if dist <= slop and dist < best_dist:
            best_dist = dist
            best_score = s
    return best_score


# ── Metric computation ────────────────────────────────────────────────────────

def _prf(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    rec  = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1   = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0.0
    return round(prec, 4), round(rec, 4), round(f1, 4)


def _auc_pr(scores: list[float], labels: list[int]) -> float:
    """
    Trapezoid-rule AUC-PR (no sklearn dependency).
    Handles ties by using micro-averaged precision at each threshold.
    """
    n_pos = sum(labels)
    if n_pos == 0 or len(set(labels)) < 2:
        return float("nan")

    # Sort by descending score
    pairs = sorted(zip(scores, labels), key=lambda x: (-x[0], -x[1]))
    tp = fp = 0
    precisions = [1.0]
    recalls    = [0.0]

    for _, label in pairs:
        if label == 1:
            tp += 1
        else:
            fp += 1
        precisions.append(tp / (tp + fp))
        recalls.append(tp / n_pos)

    # Trapezoid rule
    auc = sum(
        (recalls[i] - recalls[i - 1]) * (precisions[i] + precisions[i - 1]) / 2
        for i in range(1, len(recalls))
    )
    return round(auc, 4)


def evaluate(
    method_ids: set[str],
    truth: pd.DataFrame,
    match_slop: int,
    scores: dict[str, float] | None = None,
) -> dict:
    """
    Evaluate a predicted set against binary ground truth.
    Excludes ambiguous entries (is_true == -1).
    """
    gt = truth[truth["is_true"].isin([0, 1])]
    tp = fp = fn = tn = 0
    score_list: list[float] = []
    label_list: list[int]   = []

    for _, row in gt.iterrows():
        detected = _coord_match_in_set(row["circ_id"], method_ids, match_slop)
        lbl = int(row["is_true"])

        # Score for AUC-PR: use confidence/BSJ if available, else binary detection flag
        if scores is not None:
            sc = _score_for_id(row["circ_id"], scores, method_ids, match_slop)
        else:
            sc = 1.0 if detected else 0.0
        score_list.append(sc)
        label_list.append(lbl)

        if lbl == 1 and detected:
            tp += 1
        elif lbl == 0 and detected:
            fp += 1
        elif lbl == 1 and not detected:
            fn += 1
        else:  # lbl == 0 and not detected
            tn += 1

    prec, rec, f1 = _prf(tp, fp, fn)
    auc = _auc_pr(score_list, label_list)
    specificity = round(tn / (tn + fp), 4) if (tn + fp) > 0 else 0.0
    return {
        "n_detected": len(method_ids),
        "TP": tp, "FP": fp, "FN": fn, "TN": tn,
        "Precision": prec, "Recall": rec, "F1": f1,
        "Specificity": specificity, "AUC_PR": auc,
    }


def fp_score_binned(
    method_ids: set[str],
    summary_scores: dict[str, float],
    truth: pd.DataFrame,
    match_slop: int,
) -> list[dict]:
    """
    Bin detected circRNAs by confidence_score and label each as TP/FP.
    Returns one row per (bin, label) with count.
    Used to compare FP distributions between methods.
    """
    gt_pos = set(truth[truth["is_true"] == 1]["circ_id"])
    gt_neg = set(truth[truth["is_true"] == 0]["circ_id"])

    bins   = [0, 1.5, 2.0, 2.5, 3.0, 3.5, float("inf")]
    labels = ["<1.5", "1.5–2.0", "2.0–2.5", "2.5–3.0", "3.0–3.5", "≥3.5"]
    tp_counts = {lbl: 0 for lbl in labels}
    fp_counts = {lbl: 0 for lbl in labels}

    for circ_id, score in summary_scores.items():
        # assign bin
        bin_lbl = labels[-1]
        for i in range(len(bins) - 1):
            if bins[i] <= score < bins[i + 1]:
                bin_lbl = labels[i]
                break

        if _coord_match_in_set(circ_id, gt_pos, match_slop):
            tp_counts[bin_lbl] += 1
        elif _coord_match_in_set(circ_id, gt_neg, match_slop):
            fp_counts[bin_lbl] += 1

    rows = []
    for lbl in labels:
        tp, fp = tp_counts[lbl], fp_counts[lbl]
        fp_rate = round(fp / (tp + fp), 4) if (tp + fp) > 0 else 0.0
        rows.append({
            "score_bin": lbl,
            "TP": tp, "FP": fp, "FP_rate": fp_rate,
        })
    return rows


def stratified_f1(
    method_ids: set[str],
    truth: pd.DataFrame,
    match_slop: int,
    q1_cutoff: float,
    q3_cutoff: float,
) -> dict[str, float]:
    """F1 score stratified by BSJ count quartile in total RNA sample.

    Tier boundaries are computed from the ground-truth set (CirComPara2 / nf-core style):
      low_Q1   : bsj_total ≤ Q1  (bottom 25%)
      mid_Q2Q3 : Q1 < bsj_total ≤ Q3  (middle 50%)
      high_Q4  : bsj_total > Q3  (top 25%)
    """
    gt = truth[truth["is_true"].isin([0, 1])]
    tiers = [
        (0,          q1_cutoff, "low_Q1"),
        (q1_cutoff,  q3_cutoff, "mid_Q2Q3"),
        (q3_cutoff,  1e9,       "high_Q4"),
    ]
    result: dict[str, float] = {}
    for lo, hi, label in tiers:
        tier = gt[(gt["bsj_total"] > lo) & (gt["bsj_total"] <= hi)]
        tp = fp = fn = tn = 0
        for _, row in tier.iterrows():
            detected = _coord_match_in_set(row["circ_id"], method_ids, match_slop)
            lbl = int(row["is_true"])
            if lbl == 1 and detected:     tp += 1
            elif lbl == 0 and detected:   fp += 1
            elif lbl == 1:                fn += 1
            else:                         tn += 1
        prec, _, f1 = _prf(tp, fp, fn)
        spec = round(tn / (tn + fp), 4) if (tn + fp) > 0 else float("nan")
        result[label] = f1
        result[f"prec_{label}"] = prec
        result[f"spec_{label}"] = spec
    return result


# ── Loaders ───────────────────────────────────────────────────────────────────

def _load_bed(path: str) -> set[str]:
    ids: set[str] = set()
    with open(path) as fh:
        for line in fh:
            p = line.strip().split("\t")
            if len(p) >= 3:
                ids.add(f"{p[0]}:{p[1]}|{p[2]}")
    return ids


def _load_summary_scores(path: str) -> dict[str, float]:
    df = pd.read_csv(path, sep="\t")
    if "circ_id" in df.columns and "confidence_score" in df.columns:
        return dict(zip(df["circ_id"], df["confidence_score"]))
    return {}


def _load_dcc(path: str, min_bsj: int) -> tuple[set[str], dict[str, float]]:
    """Returns (set of circ_ids, {circ_id: count}) from DCC output.

    DCC 0.5.0 produces two files:
      CircCoordinates – annotation (Chr/Start/End/Gene/…), NO count column
      CircRNACount    – counts (Chr/Start/End/Chimeric.out.junction)
    CircCoordinates column "c" does not exist → cnt_i falls back to index 3
    (= Gene, a string) → float() fails → 0 detections.
    Fix: prefer CircRNACount (same directory) which has the actual BSJ counts.
    """
    from pathlib import Path as _Path
    count_file = _Path(path).parent / "CircRNACount"
    read_path  = str(count_file) if count_file.exists() else path

    ids: set[str]            = set()
    scores: dict[str, float] = {}
    with open(read_path) as fh:
        raw_hdr = fh.readline().strip().split("\t")
        hdr     = [h.lower().strip() for h in raw_hdr]
        chr_i   = next((i for i, h in enumerate(hdr) if h == "chr"),   0)
        start_i = next((i for i, h in enumerate(hdr) if h == "start"), 1)
        end_i   = next((i for i, h in enumerate(hdr) if h == "end"),   2)
        # Accept "c", "chimeric.out.junction", or any 4th column as count
        cnt_i   = next(
            (i for i, h in enumerate(hdr)
             if h not in ("chr", "start", "end") and i >= 3),
            3,
        )
        for line in fh:
            p = line.strip().split("\t")
            if len(p) <= max(chr_i, start_i, end_i, cnt_i):
                continue
            try:
                cnt = float(p[cnt_i])
                if cnt >= min_bsj:
                    cid = f"{p[chr_i]}:{p[start_i]}|{p[end_i]}"
                    ids.add(cid)
                    scores[cid] = cnt
            except (ValueError, IndexError):
                continue
    return ids, scores


# ── Threshold-based PR curve (multi-method) ───────────────────────────────────

def _parse_ciri2(path: str) -> dict:
    """CIRI2 output: col 1 = chr:start|end (1-based), col 5 = #junction_reads."""
    bsj = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 5:
                continue
            try:
                bsj[p[0]] = float(p[4])
            except (ValueError, IndexError):
                pass
    return bsj


def _parse_dcc_count(path: str) -> dict:
    """DCC CircRNACount: chr, start, end, count (tab-sep, 1-based coords)."""
    bsj = {}
    with open(path) as f:
        next(f)  # skip header "Chr Start End Chimeric.out.junction"
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 4:
                continue
            try:
                bsj[f"{p[0]}:{p[1]}|{p[2]}"] = float(p[3])
            except (ValueError, IndexError):
                pass
    return bsj


def _parse_ce2_binary(path: str) -> dict:
    """CIRCexplorer2 known_circ.txt (BED, 0-based start).
    Score column is unreliable (mostly 0); treat all detections as count=1.
    Start is converted to 1-based (+1) to align with CIRI2 coordinates."""
    bsj = {}
    with open(path) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) < 3 or not p[0]:
                continue
            try:
                # BED 0-based start → 1-based to match CIRI2
                cid = f"{p[0]}:{int(p[1])+1}|{p[2]}"
                bsj[cid] = 1.0
            except (ValueError, IndexError):
                pass
    return bsj


def _parse_find_circ_count(path: str) -> dict:
    """find_circ splice_sites.bed, CIRCULAR only (col 18), col 5 = junction count.
    BED 0-based start → 1-based to match CIRI2 coordinates."""
    bsj = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) < 18 or "CIRCULAR" not in p[17]:
                continue
            try:
                cid = f"{p[0]}:{int(p[1])+1}|{p[2]}"
                bsj[cid] = float(p[4])
            except (ValueError, IndexError):
                pass
    return bsj


def _build_chr_idx(bsj_dict: dict) -> dict:
    """Build {chr: [(start, end, count, cid)]} for fast coordinate lookup."""
    idx = {}
    for cid, cnt in bsj_dict.items():
        parsed = _parse_id(cid)
        if parsed:
            ch, s, e = parsed
            idx.setdefault(ch, []).append((s, e, cnt, cid))
    return idx


def _any_match(ch: str, s: int, e: int, idx: dict, thr: float, slop: int) -> bool:
    """Return True if any entry in idx[ch] passes threshold and is within slop."""
    for ds, de, dc, _ in idx.get(ch, []):
        if dc >= thr and max(abs(ds - s), abs(de - e)) <= slop:
            return True
    return False


def _evaluate_consensus(consensus: set, gt: "pd.DataFrame", slop: int) -> tuple:
    tp = fp = fn = tn = 0
    for _, row in gt.iterrows():
        detected = _coord_match_in_set(row["circ_id"], consensus, slop)
        lbl = int(row["is_true"])
        if lbl == 1 and detected:       tp += 1
        elif lbl == 0 and detected:     fp += 1
        elif lbl == 1 and not detected: fn += 1
        else:                           tn += 1
    return tp, fp, fn, tn


def _auc_from_pr_points(pr_points: list) -> float:
    pr_sorted = sorted(set(pr_points))
    if not pr_sorted or pr_sorted[0][0] > 0:
        pr_sorted.insert(0, (0.0, 1.0))
    return round(sum(
        (pr_sorted[i][0] - pr_sorted[i-1][0])
        * (pr_sorted[i][1] + pr_sorted[i-1][1]) / 2
        for i in range(1, len(pr_sorted))
    ), 4)


def _multi_method_pr_curves(
    ciri2_path: str,
    dcc_count_path: str,
    ce2_path: str,
    find_circ_path: str,
    truth: "pd.DataFrame",
    thresholds: list,
    slop_our: int = 10,
    slop_cp2: int = 10,
    slop_nfc: int = 1,  # 1 instead of 0 to absorb BED→1-based +1 offset residuals
) -> dict:
    """
    Sweep min_bsj threshold across 3 methods simultaneously.
    All methods seeded from CIRI2 (primary tool). Secondary support checked:
      Our_adaptive:       CIRI2 AND DCC
      CirComPara2_4tools: CIRI2 AND (DCC OR CE2 OR find_circ)
      nfcore_3tools:      CIRI2 AND (CE2 OR find_circ)  [CE2 binary, slop=1]

    CIRCexplorer2 is treated as binary (all detections, no count threshold)
    because its score column is unreliable in this dataset (97% zeros).

    Returns: {method_name: (auc_pr, [row_dicts])}
    """
    print("[pr_curve] Parsing tool files...", file=sys.stderr)
    ciri2_bsj  = _parse_ciri2(ciri2_path) if ciri2_path else {}
    dcc_bsj    = _parse_dcc_count(dcc_count_path) if dcc_count_path else {}
    ce2_bsj    = _parse_ce2_binary(ce2_path) if ce2_path else {}    # all count=1
    fc_bsj     = _parse_find_circ_count(find_circ_path) if find_circ_path else {}

    dcc_idx = _build_chr_idx(dcc_bsj)
    ce2_idx = _build_chr_idx(ce2_bsj)
    fc_idx  = _build_chr_idx(fc_bsj)

    print(f"[pr_curve] CIRI2={len(ciri2_bsj)}, DCC={len(dcc_bsj)}, "
          f"CE2={len(ce2_bsj)}, find_circ={len(fc_bsj)}", file=sys.stderr)

    gt = truth[truth["is_true"].isin([0, 1])]

    method_names = ["Our_adaptive", "CirComPara2_4tools", "nfcore_3tools"]
    pr_pts  = {m: [] for m in method_names}
    all_rows = {m: [] for m in method_names}

    for thr in thresholds:
        c2_ids = {cid for cid, bsj in ciri2_bsj.items() if bsj >= thr}

        our_con = set()
        cp2_con = set()
        nfc_con = set()

        for cid in c2_ids:
            parsed = _parse_id(cid)
            if not parsed:
                continue
            ch, s, e = parsed

            dcc_ok  = _any_match(ch, s, e, dcc_idx, thr,  slop_our)
            dcc_ok2 = _any_match(ch, s, e, dcc_idx, thr,  slop_cp2)
            ce2_ok  = _any_match(ch, s, e, ce2_idx, 1.0,  slop_cp2)
            ce2_nfc = _any_match(ch, s, e, ce2_idx, 1.0,  slop_nfc)
            fc_ok   = _any_match(ch, s, e, fc_idx,  thr,  slop_cp2)
            fc_nfc  = _any_match(ch, s, e, fc_idx,  thr,  slop_nfc)

            if dcc_ok:
                our_con.add(cid)
            if dcc_ok2 or ce2_ok or fc_ok:
                cp2_con.add(cid)
            if ce2_nfc or fc_nfc:
                nfc_con.add(cid)

        for name, con, slop in [
            ("Our_adaptive",       our_con, slop_our),
            ("CirComPara2_4tools", cp2_con, slop_cp2),
            ("nfcore_3tools",      nfc_con, slop_nfc),
        ]:
            tp, fp, fn, _ = _evaluate_consensus(con, gt, slop)
            prec, rec, f1 = _prf(tp, fp, fn)
            all_rows[name].append({
                "threshold": thr, "n_detected": len(con),
                "TP": tp, "FP": fp,
                "Precision": prec, "Recall": rec, "F1": f1,
            })
            pr_pts[name].append((rec, prec))
            print(f"  [{name:22s}] bsj>={thr:3d}: n={len(con):5d}  "
                  f"Prec={prec:.3f}  Rec={rec:.3f}  F1={f1:.3f}", file=sys.stderr)

    return {m: (_auc_from_pr_points(pr_pts[m]), all_rows[m]) for m in method_names}


# ── Legacy single-method wrapper (kept for backward compat) ───────────────────

def _threshold_pr_curve(
    ciri2_path: str,
    dcc_count_path: str,
    truth: "pd.DataFrame",
    thresholds: list,
    slop: int = 10,
) -> tuple:
    """Legacy wrapper: Our_adaptive only. Use _multi_method_pr_curves for new code."""
    result = _multi_method_pr_curves(
        ciri2_path, dcc_count_path, "", "", truth, thresholds, slop_our=slop,
    )
    return result["Our_adaptive"]


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Evaluate circRNA detection accuracy against RNase R ground truth"
    )
    parser.add_argument("--ground-truth",          required=True)
    parser.add_argument("--our-bed",               required=True,
                        help="Our consensus filter BED (consensus_filter.py output)")
    parser.add_argument("--our-summary",           required=True,
                        help="Consensus summary TSV (contains confidence_score)")
    parser.add_argument("--our-no-qc-bed",         default=None,
                        help="Our method WITHOUT pseudo-circ QC (max_junction_ratio=99)")
    parser.add_argument("--our-no-qc-summary",     default=None,
                        help="Summary TSV for no-QC ablation")
    parser.add_argument("--circompara2-bed",        default=None,
                        help="CirComPara2 sim BED (slop=10, no BSJ/FSJ QC; optional ablation)")
    parser.add_argument("--circompara2-summary",    default=None,
                        help="CirComPara2 sim summary TSV")
    parser.add_argument("--circompara2-4tools-bed",     default=None,
                        help="CirComPara2 full 5-tool BED (CIRI2+CIRIquant+DCC+CIRCexplorer2+find_circ, slop=10)")
    parser.add_argument("--circompara2-4tools-summary", default=None,
                        help="CirComPara2 5-tool summary TSV")
    parser.add_argument("--nfcore-bed",            required=True,
                        help="nf-core 3-tool sim BED (CIRIquant+CIRCexplorer2+find_circ, slop=0)")
    parser.add_argument("--nfcore-summary",        default=None,
                        help="nf-core 3-tool sim summary TSV (contains confidence_score)")
    parser.add_argument("--output-summary",        required=True)
    parser.add_argument("--output-stratified",     required=True)
    parser.add_argument("--output-fp-comparison",  default=None,
                        dest="output_fp_comparison",
                        help="TSV: FP score distribution comparison (Our vs CirComPara2)")
    parser.add_argument("--output-pr-curve", default=None, dest="output_pr_curve",
                        help="TSV: threshold-based PR curve for all 3 methods (honest AUC-PR)")
    parser.add_argument("--ciri2-file", default=None, dest="ciri2_file",
                        help="CIRI2 output file (*.ciri2) for threshold sweep")
    parser.add_argument("--dcc-count-file", default=None, dest="dcc_count_file",
                        help="DCC CircRNACount file for threshold sweep")
    parser.add_argument("--circexplorer2-file", default=None, dest="circexplorer2_file",
                        help="CIRCexplorer2 known_circ.txt for threshold sweep (treated as binary)")
    parser.add_argument("--find-circ-file", default=None, dest="find_circ_file",
                        help="find_circ splice_sites.bed for threshold sweep (CIRCULAR only)")
    parser.add_argument("--slop",    type=int, default=10)
    parser.add_argument("--min-bsj", type=int, default=2, dest="min_bsj")
    args = parser.parse_args()

    truth = pd.read_csv(args.ground_truth, sep="\t")
    n_tp  = int((truth["is_true"] == 1).sum())
    n_tn  = int((truth["is_true"] == 0).sum())
    print(f"[accuracy] Ground truth: TP={n_tp}, TN={n_tn} "
          f"(ambiguous excluded)", file=sys.stderr)

    # Quartile boundaries for stratified F1 (CirComPara2 / nf-core style)
    _gt_bsj = truth.loc[truth["is_true"].isin([0, 1]), "bsj_total"].dropna()
    q1_cutoff = float(np.percentile(_gt_bsj, 25))
    q3_cutoff = float(np.percentile(_gt_bsj, 75))
    print(f"[accuracy] BSJ count quartiles: Q1={q1_cutoff:.1f}, Q3={q3_cutoff:.1f}",
          file=sys.stderr)

    # ── Load predictions ──────────────────────────────────────────────────────
    our_ids              = _load_bed(args.our_bed)
    our_scores           = _load_summary_scores(args.our_summary)
    our_no_qc_ids        = _load_bed(args.our_no_qc_bed) if args.our_no_qc_bed else None
    our_no_qc_scores     = _load_summary_scores(args.our_no_qc_summary) if args.our_no_qc_summary else None
    circompara2_ids    = _load_bed(args.circompara2_bed) if args.circompara2_bed else None
    circompara2_scores = _load_summary_scores(args.circompara2_summary) if args.circompara2_summary else None
    cp2_4t_ids    = _load_bed(args.circompara2_4tools_bed) if args.circompara2_4tools_bed else None
    cp2_4t_scores = _load_summary_scores(args.circompara2_4tools_summary) if args.circompara2_4tools_summary else None
    nfcore_ids     = _load_bed(args.nfcore_bed)
    nfcore_scores  = (
        _load_summary_scores(args.nfcore_summary)
        if args.nfcore_summary else None
    )

    det_msg = f"[accuracy] Detected: ours={len(our_ids)}"
    if our_no_qc_ids is not None:
        det_msg += f", ours_no_qc={len(our_no_qc_ids)}"
    if circompara2_ids is not None:
        det_msg += f", circompara2_sim={len(circompara2_ids)}"
    if cp2_4t_ids is not None:
        det_msg += f", circompara2_4tools={len(cp2_4t_ids)}"
    det_msg += f", nfcore_3tools={len(nfcore_ids)}"
    print(det_msg, file=sys.stderr)

    # ── Evaluate each method ──────────────────────────────────────────────────
    methods = [("Our_adaptive", our_ids, args.slop, our_scores)]
    if our_no_qc_ids is not None:
        methods.append(("Our_no_QC", our_no_qc_ids, args.slop, our_no_qc_scores))
    if circompara2_ids is not None:
        methods.append(("CirComPara2_sim", circompara2_ids, args.slop, circompara2_scores))
    if cp2_4t_ids is not None:
        methods.append(("CirComPara2_4tools", cp2_4t_ids, args.slop, cp2_4t_scores))
    methods.append(("nfcore_3tools", nfcore_ids, 0, nfcore_scores))

    summary_rows = []
    strat_rows   = []
    for name, ids, slop, scores in methods:
        m = evaluate(ids, truth, slop, scores)
        summary_rows.append({"Method": name, **m})

        s = stratified_f1(ids, truth, slop, q1_cutoff, q3_cutoff)
        strat_rows.append({"Method": name, **s,
                           "q1_cutoff": q1_cutoff, "q3_cutoff": q3_cutoff})
        print(
            f"[accuracy] {name:22s}  "
            f"Prec={m['Precision']:.3f}  Rec={m['Recall']:.3f}  "
            f"F1={m['F1']:.3f}  Spec={m['Specificity']:.3f}  "
            f"TN={m['TN']}  AUC-PR={m['AUC_PR']:.3f}",
            file=sys.stderr,
        )

    # ── Threshold-based PR curve + honest AUC-PR (all 3 methods) ─────────────
    if args.ciri2_file and args.dcc_count_file:
        PR_THRESHOLDS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50]
        print("[accuracy] Computing threshold-based PR curves (all 3 methods)...",
              file=sys.stderr)
        multi_results = _multi_method_pr_curves(
            args.ciri2_file,
            args.dcc_count_file,
            args.circexplorer2_file or "",
            args.find_circ_file or "",
            truth,
            PR_THRESHOLDS,
            slop_our=args.slop,
            slop_cp2=args.slop,
            slop_nfc=1,  # nfcore: slop=0 nominal, +1 for BED→1-based offset
        )
        # Override AUC_PR in summary rows
        method_map = {
            "Our_adaptive":       "Our_adaptive",
            "CirComPara2_4tools": "CirComPara2_4tools",
            "nfcore_3tools":      "nfcore_3tools",
        }
        for row in summary_rows:
            m = row["Method"]
            if m in method_map and method_map[m] in multi_results:
                auc, _ = multi_results[method_map[m]]
                old_auc = row["AUC_PR"]
                row["AUC_PR"] = auc
                row["AUC_PR_note"] = "threshold-sweep"
                print(f"[accuracy] {m}: AUC-PR {old_auc:.4f} → {auc:.4f} (threshold-sweep)",
                      file=sys.stderr)
            else:
                row["AUC_PR_note"] = "binary-detection"

        if args.output_pr_curve:
            # Long format: threshold, method, n_detected, TP, FP, Precision, Recall, F1
            pr_long = []
            for m, (_, rows) in multi_results.items():
                for r in rows:
                    pr_long.append({"method": m, **r})
            pd.DataFrame(pr_long).to_csv(args.output_pr_curve, sep="\t", index=False)
            print(f"[accuracy] PR curve (3 methods) → {args.output_pr_curve}", file=sys.stderr)

    # ── Write outputs ─────────────────────────────────────────────────────────
    out_dir = Path(args.output_summary).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(summary_rows).to_csv(args.output_summary,    sep="\t", index=False)
    pd.DataFrame(strat_rows).to_csv(  args.output_stratified, sep="\t", index=False)
    print(f"[accuracy] Summary    → {args.output_summary}",    file=sys.stderr)
    print(f"[accuracy] Stratified → {args.output_stratified}", file=sys.stderr)

    # ── FP score distribution comparison (Our vs CirComPara2) ─────────────────
    cp2_ref_ids    = circompara2_ids    if circompara2_ids    is not None else cp2_4t_ids
    cp2_ref_scores = circompara2_scores if circompara2_scores is not None else cp2_4t_scores
    if args.output_fp_comparison and cp2_ref_ids is not None:
        our_bins = fp_score_binned(our_ids, our_scores, truth, args.slop)
        cp2_bins = fp_score_binned(cp2_ref_ids, cp2_ref_scores, truth, args.slop)
        fp_rows = []
        for o, c in zip(our_bins, cp2_bins):
            fp_rows.append({
                "score_bin":          o["score_bin"],
                "Our_TP":             o["TP"],
                "Our_FP":             o["FP"],
                "Our_FP_rate":        o["FP_rate"],
                "CirComPara2_TP":     c["TP"],
                "CirComPara2_FP":     c["FP"],
                "CirComPara2_FP_rate":c["FP_rate"],
            })
        pd.DataFrame(fp_rows).to_csv(args.output_fp_comparison, sep="\t", index=False)
        print(f"[accuracy] FP comparison → {args.output_fp_comparison}", file=sys.stderr)


if __name__ == "__main__":
    main()
