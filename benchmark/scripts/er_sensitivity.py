"""ER-threshold and replicate sensitivity of the RNase R ground truth (benchmark).

For each variant (full = both RNase R replicates averaged, as in rnaser_ground_truth.py;
rep1_only = replicate 1 only) and each threshold pair, rebuilds the TP/TN/ambiguous labels on
the same total-RNA universe, then scores the three detection strategies with
accuracy_benchmark.evaluate (same slop and exclusion of ambiguous entries).

The labelling logic is the same as rnaser_ground_truth.main(); the reproduction of the
published ground truth (full, 1.5 / 0.5) is checked row by row before any other row is written.

Threshold grid: positive threshold p in {1.2, 1.5, 2, 3, 5} with negative threshold 1/p, plus the
asymmetric reference pair 1.5 / 0.5 used by the paper. Note 1/1.5 = 0.667 != 0.5: the reference
pair is not a reciprocal pair, and the output marks it.

usage:
  python3 er_sensitivity.py --total <total.gtf> --rep1 <rep1.gtf> --rep2 <rep2.gtf> \
      --published-gt <ground_truth.tsv> --beds <circDEX.bed> <CirComPara2.bed> <nfcore.bed> \
      --output <er_sensitivity.tsv>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import accuracy_benchmark as ab  # noqa: E402
import rnaser_ground_truth as rgt  # noqa: E402

THRESHOLDS = [(1.2, round(1 / 1.2, 4)), (1.5, round(1 / 1.5, 4)), (2.0, 0.5),
              (3.0, round(1 / 3, 4)), (5.0, 0.2), (1.5, 0.5)]
STRATEGIES = ["circDEX", "CirComPara2_4tools", "nfcore_3tools"]
SLOP = 10


def build_gt(total_recs: dict, rnaser_counts: list, er_pos: float, er_neg: float, slop: int = SLOP) -> pd.DataFrame:
    total_counts = rgt._bsj_counts(total_recs)
    total_index = rgt._build_index(total_counts)
    rnaser_indices = [rgt._build_index(c) for c in rnaser_counts]
    rows = []
    for circ_id in sorted(total_recs):
        bsj_total = rgt._coord_match(circ_id, total_index, slop)
        vals = [rgt._coord_match(circ_id, idx, slop) for idx in rnaser_indices]
        bsj_rnaser_mean = sum(vals) / len(vals)
        if bsj_total < 0.001:
            er = 100.0 if bsj_rnaser_mean > 0.001 else 1.0
        else:
            er = bsj_rnaser_mean / bsj_total
        if er >= er_pos:
            is_true = 1
        elif er <= er_neg:
            is_true = 0
        else:
            is_true = -1
        m = total_recs[circ_id]
        rows.append({"circ_id": circ_id, "chr": m["chr"], "start": m["start"], "end": m["end"],
                     "bsj_total": round(bsj_total, 4), "bsj_rnaser_mean": round(bsj_rnaser_mean, 4),
                     "enrichment_ratio": round(er, 4), "is_true": is_true})
    return pd.DataFrame(rows)


def validate_against_published(gt: pd.DataFrame, published_path: str) -> None:
    pub = pd.read_csv(published_path, sep="\t")
    a = gt.set_index("circ_id")[["is_true", "enrichment_ratio"]]
    b = pub.set_index("circ_id")[["is_true", "enrichment_ratio"]]
    if len(a) != len(b) or set(a.index) != set(b.index):
        raise SystemExit(f"[er] VALIDATION FAILED: row sets differ ({len(a)} vs {len(b)})")
    j = a.join(b, lsuffix="_new", rsuffix="_pub")
    n_lbl = int((j.is_true_new != j.is_true_pub).sum())
    d_er = float((j.enrichment_ratio_new - j.enrichment_ratio_pub).abs().max())
    if n_lbl or d_er > 1e-4:
        raise SystemExit(f"[er] VALIDATION FAILED: {n_lbl} label differences, max |dER|={d_er:.3g}")
    print(f"[er] validation OK: {len(j)} rows, 0 label differences, max |dER|={d_er:.2g}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--total", required=True)
    ap.add_argument("--rep1", required=True)
    ap.add_argument("--rep2", required=True)
    ap.add_argument("--published-gt", required=True, dest="published_gt")
    ap.add_argument("--beds", nargs=3, required=True, help="circDEX, CirComPara2_4tools, nfcore_3tools BED")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    total_recs = rgt._parse_gtf(args.total)
    rep1 = rgt._parse_gtf(args.rep1)
    rep2 = rgt._parse_gtf(args.rep2)
    rep1_c = rgt._bsj_counts(rep1)
    rep2_c = rgt._bsj_counts(rep2)
    detections = {name: ab._load_bed(p) for name, p in zip(STRATEGIES, args.beds)}

    variants = {"full": [rep1_c, rep2_c], "rep1_only": [rep1_c]}
    # full, 1.5/0.5: exact reproduction of the published ground truth
    gt_full = build_gt(total_recs, variants["full"], 1.5, 0.5)
    validate_against_published(gt_full, args.published_gt)

    rows = []
    for vname, rcounts in variants.items():
        for pos, neg in THRESHOLDS:
            gt = gt_full if (vname == "full" and (pos, neg) == (1.5, 0.5)) else build_gt(total_recs, rcounts, pos, neg)
            n_pos = int((gt.is_true == 1).sum())
            n_neg = int((gt.is_true == 0).sum())
            n_amb = int((gt.is_true == -1).sum())
            row = {"variant": vname, "er_pos": pos, "er_neg": neg,
                   "reciprocal_pair": "yes" if abs(pos * neg - 1) < 1e-3 else "no (reference)",
                   "n_reference": len(gt), "n_positive": n_pos, "n_negative": n_neg, "n_excluded": n_amb}
            for s in STRATEGIES:
                ev = ab.evaluate(detections[s], gt, SLOP)
                row[f"{s}_precision"] = ev["Precision"]
                row[f"{s}_recall"] = ev["Recall"]
                row[f"{s}_specificity"] = ev["Specificity"]
            rows.append(row)
            print(f"[er] {vname} {pos}/{neg}: pos={n_pos} neg={n_neg} amb={n_amb}", file=sys.stderr)

    df = pd.DataFrame(rows)
    for metric in ["precision", "recall", "specificity"]:
        cols = [f"{s}_{metric}" for s in STRATEGIES]
        ranks = df[cols].rank(axis=1, ascending=False, method="min").astype(int)
        for s in STRATEGIES:
            df[f"rank_{metric}_{s}"] = ranks[f"{s}_{metric}"]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[er] written → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
