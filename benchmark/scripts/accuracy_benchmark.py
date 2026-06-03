"""
accuracy_benchmark.py – 偵測準確率比較（對照 RNase R ground truth）

Five detection strategies evaluated:
  1. Our adaptive consensus filter
     CIRIquant + DCC, slop=10 bp, BSJ/FSJ pseudo-circ QC
  2. CirComPara2 simulation (Gaffo et al. 2022)
     CIRIquant + DCC, slop=10 bp, NO pseudo-circ QC
  3. nf-core/circrna simulation (Digby-Bell et al. 2023)
     CIRIquant + CIRCexplorer2 + find_circ, slop=0 (exact coords), min_tools=2
  4. circRNA-sponging simulation (Hansen et al. 2023)
     DCC-only, min_bsj=2
  5. CLEAR simulation (custom 2020)
     CIRIquant-only, slop=0 (exact coords), no consensus step

Metrics: Precision, Recall, F1, AUC-PR
Stratification: low BSJ (1–4), mid (5–19), high (≥20 RPM in total RNA)

Outputs:
  --output-summary    results/benchmark/accuracy_summary.tsv
  --output-stratified results/benchmark/stratified_f1.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

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
) -> dict[str, float]:
    """F1 score stratified by BSJ RPM tier in total RNA sample."""
    gt  = truth[truth["is_true"].isin([0, 1])]
    tiers = [
        (1.0,  4.99,  "low_1-4"),
        (5.0,  19.99, "mid_5-19"),
        (20.0, 1e9,   "high_ge20"),
    ]
    result: dict[str, float] = {}
    for lo, hi, label in tiers:
        tier = gt[(gt["bsj_total"] >= lo) & (gt["bsj_total"] <= hi)]
        tp = fp = fn = 0
        for _, row in tier.iterrows():
            detected = _coord_match_in_set(row["circ_id"], method_ids, match_slop)
            lbl = int(row["is_true"])
            if lbl == 1 and detected:   tp += 1
            elif lbl == 0 and detected: fp += 1
            elif lbl == 1:              fn += 1
        _, _, f1 = _prf(tp, fp, fn)
        result[label] = f1
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
    parser.add_argument("--circompara2-bed",        required=True,
                        help="CirComPara2 sim BED (slop=10, no BSJ/FSJ QC)")
    parser.add_argument("--circompara2-summary",    required=True,
                        help="CirComPara2 sim summary TSV")
    parser.add_argument("--nfcore-bed",            required=True,
                        help="nf-core 3-tool sim BED (CIRIquant+CIRCexplorer2+find_circ, slop=0)")
    parser.add_argument("--nfcore-summary",        default=None,
                        help="nf-core 3-tool sim summary TSV (contains confidence_score)")
    parser.add_argument("--dcc-coords",            required=True,
                        help="DCC CircCoordinates (circRNA-sponging sim)")
    parser.add_argument("--clear-bed",             required=True,
                        help="CLEAR sim BED (CIRIquant-only, slop=0)")
    parser.add_argument("--output-summary",        required=True)
    parser.add_argument("--output-stratified",     required=True)
    parser.add_argument("--output-fp-comparison",  default=None,
                        dest="output_fp_comparison",
                        help="TSV: FP score distribution comparison (Our vs CirComPara2)")
    parser.add_argument("--slop",    type=int, default=10)
    parser.add_argument("--min-bsj", type=int, default=2, dest="min_bsj")
    args = parser.parse_args()

    truth = pd.read_csv(args.ground_truth, sep="\t")
    n_tp  = int((truth["is_true"] == 1).sum())
    n_tn  = int((truth["is_true"] == 0).sum())
    print(f"[accuracy] Ground truth: TP={n_tp}, TN={n_tn} "
          f"(ambiguous excluded)", file=sys.stderr)

    # ── Load predictions ──────────────────────────────────────────────────────
    our_ids              = _load_bed(args.our_bed)
    our_scores           = _load_summary_scores(args.our_summary)
    circompara2_ids      = _load_bed(args.circompara2_bed)
    circompara2_scores   = _load_summary_scores(args.circompara2_summary)
    nfcore_ids           = _load_bed(args.nfcore_bed)
    nfcore_scores        = (
        _load_summary_scores(args.nfcore_summary)
        if args.nfcore_summary else None
    )
    sponge_ids, sponge_scores = _load_dcc(args.dcc_coords, args.min_bsj)
    clear_ids            = _load_bed(args.clear_bed)

    print(
        f"[accuracy] Detected: ours={len(our_ids)}, "
        f"circompara2_sim={len(circompara2_ids)}, "
        f"nfcore_3tools={len(nfcore_ids)}, "
        f"sponging_sim={len(sponge_ids)}, "
        f"clear_sim={len(clear_ids)}",
        file=sys.stderr,
    )

    # ── Evaluate each method ──────────────────────────────────────────────────
    # slop=0 for exact-coordinate methods (nfcore, clear): evaluating at their native resolution.
    # nfcore_3tools uses confidence_score from consensus_filter.py output (no longer binary).
    methods = [
        ("Our_adaptive",    our_ids,         args.slop, our_scores),
        ("CirComPara2_sim", circompara2_ids,  args.slop, circompara2_scores),
        ("nfcore_3tools",   nfcore_ids,       0,         nfcore_scores),
        ("sponging_DCC",    sponge_ids,       args.slop, sponge_scores),
        ("CLEAR_sim",       clear_ids,        0,         None),
    ]

    summary_rows = []
    strat_rows   = []
    for name, ids, slop, scores in methods:
        m = evaluate(ids, truth, slop, scores)
        summary_rows.append({"Method": name, **m})

        s = stratified_f1(ids, truth, slop)
        strat_rows.append({"Method": name, **s})
        print(
            f"[accuracy] {name:22s}  "
            f"Prec={m['Precision']:.3f}  Rec={m['Recall']:.3f}  "
            f"F1={m['F1']:.3f}  Spec={m['Specificity']:.3f}  "
            f"TN={m['TN']}  AUC-PR={m['AUC_PR']:.3f}",
            file=sys.stderr,
        )

    # ── Write outputs ─────────────────────────────────────────────────────────
    out_dir = Path(args.output_summary).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(summary_rows).to_csv(args.output_summary,    sep="\t", index=False)
    pd.DataFrame(strat_rows).to_csv(  args.output_stratified, sep="\t", index=False)
    print(f"[accuracy] Summary    → {args.output_summary}",    file=sys.stderr)
    print(f"[accuracy] Stratified → {args.output_stratified}", file=sys.stderr)

    # ── FP score distribution comparison (Our vs CirComPara2) ─────────────────
    if args.output_fp_comparison:
        our_bins = fp_score_binned(our_ids, our_scores, truth, args.slop)
        cp2_bins = fp_score_binned(circompara2_ids, circompara2_scores, truth, args.slop)
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
