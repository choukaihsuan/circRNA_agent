"""
rnaser_ground_truth.py – Build a circRNA ground truth set from RNase R enrichment.

RNase R是一種外切核酸酶，能降解線性RNA但保留環狀RNA。
Enrichment Ratio (ER) = RNase R 處理後的 BSJ RPM / Total RNA 的 BSJ RPM
  ER ≥ er_tp → True Positive（環狀 RNA，被 RNase R 富集）
  ER ≤ er_tn → True Negative（線性偽陽性，被 RNase R 降解）
  中間值  → Ambiguous（排除）

Output columns:
  circ_id, chr, start, end, bsj_total, bsj_rnaser_mean,
  enrichment_ratio, is_true
  is_true: 1=TP, 0=TN, -1=ambiguous

Usage:
  python rnaser_ground_truth.py \
      --total  results/circRNA/SRR444655/SRR444655.gtf \
      --rnaser results/circRNA/SRR444974/SRR444974.gtf \
               results/circRNA/SRR445016/SRR445016.gtf \
      --output results/benchmark/rnaser/ground_truth.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# ── GTF parser ────────────────────────────────────────────────────────────────

def _parse_gtf(path: str) -> dict[str, dict]:
    """
    Parse CIRIquant GTF.
    Returns { circ_id : {chr, start, end, bsj, fsj} }.
    """
    records: dict[str, dict] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "circRNA":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end   = int(parts[4])
            attrs = parts[8]
            bsj_m = re.search(r'BSJ\s+([\d.]+)', attrs)
            fsj_m = re.search(r'FSJ\s+([\d.]+)', attrs)
            if not bsj_m:
                continue
            bsj = float(bsj_m.group(1))
            fsj = float(fsj_m.group(1)) if fsj_m else 0.0
            circ_id = f"{chrom}:{start}|{end}"
            records[circ_id] = {
                "chr": chrom, "start": start, "end": end,
                "bsj": bsj, "fsj": fsj,
            }
    return records


def _total_bsj(records: dict) -> float:
    return max(sum(v["bsj"] for v in records.values()), 1.0)


def _bsj_rpm(records: dict) -> dict[str, float]:
    """Convert BSJ raw counts to RPM (per million total BSJ reads)."""
    total = _total_bsj(records)
    return {k: v["bsj"] / total * 1_000_000 for k, v in records.items()}


# ── Fuzzy coordinate matching ─────────────────────────────────────────────────

def _coord_match(circ_id: str, pool: dict[str, float], slop: int) -> float:
    """
    Find the best-matching entry in pool for circ_id within slop bp.
    Returns the matched RPM value, or 0.0 if no match.
    """
    m = re.match(r'^(.+):(\d+)\|(\d+)$', circ_id)
    if not m:
        return 0.0
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))

    best_val, best_dist = 0.0, float("inf")
    for key, val in pool.items():
        km = re.match(r'^(.+):(\d+)\|(\d+)$', key)
        if not km or km.group(1) != chrom:
            continue
        dist = max(abs(int(km.group(2)) - start), abs(int(km.group(3)) - end))
        if dist <= slop and dist < best_dist:
            best_dist = dist
            best_val = val
    return best_val


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build RNase R enrichment ground truth for circRNA benchmark"
    )
    parser.add_argument(
        "--total", required=True,
        help="CIRIquant GTF for total RNA sample (e.g., SRR444655)"
    )
    parser.add_argument(
        "--rnaser", required=True, nargs="+",
        help="CIRIquant GTF(s) for RNase R-treated samples (≥1)"
    )
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument(
        "--slop", type=int, default=10,
        help="Coordinate tolerance (bp) for cross-sample matching (default: 10)"
    )
    parser.add_argument(
        "--er-tp", type=float, default=1.5, dest="er_tp",
        help="Enrichment ratio threshold for true positives (default: 1.5)"
    )
    parser.add_argument(
        "--er-tn", type=float, default=0.5, dest="er_tn",
        help="Enrichment ratio threshold for true negatives (default: 0.5)"
    )
    args = parser.parse_args()

    # ── Parse input GTFs ──────────────────────────────────────────────────────
    total_recs = _parse_gtf(args.total)
    total_rpm  = _bsj_rpm(total_recs)
    print(f"[ground_truth] Total RNA ({Path(args.total).stem}): "
          f"{len(total_recs)} circRNAs", file=sys.stderr)

    rnaser_rpms: list[dict[str, float]] = []
    for path in args.rnaser:
        recs = _parse_gtf(path)
        rnaser_rpms.append(_bsj_rpm(recs))
        print(f"[ground_truth] RNase R  ({Path(path).stem}): "
              f"{len(recs)} circRNAs", file=sys.stderr)

    # ── Union of all detected circRNAs ────────────────────────────────────────
    all_ids: set[str] = set(total_recs.keys())
    for rpm_map in rnaser_rpms:
        all_ids.update(rpm_map.keys())

    # ── Compute enrichment ratios ─────────────────────────────────────────────
    rows = []
    for circ_id in sorted(all_ids):
        bsj_total  = _coord_match(circ_id, total_rpm, args.slop)
        bsj_rnaser_vals = [
            _coord_match(circ_id, rpm_map, args.slop) for rpm_map in rnaser_rpms
        ]
        bsj_rnaser_mean = sum(bsj_rnaser_vals) / len(bsj_rnaser_vals)

        # Enrichment ratio with edge-case handling
        if bsj_total < 0.001:
            if bsj_rnaser_mean > 0.001:
                er = 100.0   # only in RNase R → strong TP signal
            else:
                er = 1.0     # not detected in either → ambiguous
        else:
            er = bsj_rnaser_mean / bsj_total

        # Label: 1=TP, 0=TN, -1=ambiguous
        if er >= args.er_tp:
            is_true = 1
        elif er <= args.er_tn:
            is_true = 0
        else:
            is_true = -1

        # Coordinates come from whichever GTF detected this circRNA
        meta = total_recs.get(circ_id)
        if meta is None:
            # Fall back to RNase R GTF coordinates
            m = re.match(r'^(.+):(\d+)\|(\d+)$', circ_id)
            meta = {"chr": m.group(1), "start": int(m.group(2)), "end": int(m.group(3))} \
                   if m else {"chr": ".", "start": 0, "end": 0}

        rows.append({
            "circ_id":          circ_id,
            "chr":              meta["chr"],
            "start":            meta["start"],
            "end":              meta["end"],
            "bsj_total":        round(bsj_total, 4),
            "bsj_rnaser_mean":  round(bsj_rnaser_mean, 4),
            "enrichment_ratio": round(er, 4),
            "is_true":          is_true,
        })

    df = pd.DataFrame(rows)

    n_tp  = int((df["is_true"] == 1).sum())
    n_tn  = int((df["is_true"] == 0).sum())
    n_amb = int((df["is_true"] == -1).sum())
    print(
        f"[ground_truth] ER threshold: TP ≥ {args.er_tp}, TN ≤ {args.er_tn}\n"
        f"[ground_truth] TP={n_tp}  TN={n_tn}  ambiguous={n_amb}  "
        f"total={len(df)}",
        file=sys.stderr,
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[ground_truth] Written → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
