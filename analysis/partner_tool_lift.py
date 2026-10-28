#!/usr/bin/env python3
"""
Item 15 Stage 0 (reviewer Q2 follow-up): for each cell line, stratify
CIRIquant's own ground-truth-evaluable calls by HOW MANY of the three
candidate partner tools (DCC, CIRCexplorer2, find_circ) also call that
locus (0-3), and report per-stratum locus count, TP rate, and lift
(stratum TP rate / overall TP rate) -- not precision, since these
partners are not ground truth and a precision framing would imply they
are.

Read-only throughout: reuses consensus_filter.py's own parser functions
(imported, not reimplemented) and each cell line's existing
ground_truth.tsv / CIRIquant GTF / DCC CircCoordinates / CIRCexplorer2
known_circ.txt / find_circ splice_sites.bed. No detection tool rerun, no
Snakemake, no new consensus computation.

Per-tool slop, recorded per CLI dataset rather than as one global
constant (per task -- this is the actual root cause behind a separate,
earlier B1 finding: nf-core's slop=0 is its own published strategy, not a
bug to be "fixed" by unifying to circDEX's 10bp):
  DCC           -> 10bp  (circDEX's own project-wide consensus convention;
                           DCC never appears anywhere in this project's
                           benchmarks at exact-match)
  CIRCexplorer2 -> 0bp   (nf-core/circrna's own published workflow uses
                           exact coordinate matching for this tool)
  find_circ     -> 0bp   (same reasoning -- one of nf-core/circrna's three
                           tools, benchmarked at nf-core's own convention)
  CIRIquant-vs-ground-truth classification -> 10bp (circDEX's own,
                           project-wide convention; this is about
                           classifying CIRIquant's OWN calls, not about
                           partner corroboration, so circDEX's own slop
                           applies here regardless of which partner is
                           being discussed)
"""
import bisect
import re
import sys
import pandas as pd

sys.path.insert(0, "/home/choukaihsuan/circRNA_agent/scripts")
import consensus_filter as cf

MIN_BSJ = 2
SLOP_GT = 10          # CIRIquant vs ground truth (circDEX's own convention)
PARTNER_SLOP = {"DCC": 10, "CIRCexplorer2": 0, "find_circ": 0}

DATASETS = {
    "Hs68": dict(
        ciriquant_gtf="/home/choukaihsuan/GSE55872_results/circRNA/SRR444655/SRR444655.gtf",
        dcc="/home/choukaihsuan/GSE55872_results/circRNA/SRR444655/DCC/CircCoordinates",
        ce2="/home/choukaihsuan/GSE55872_results/circRNA/SRR444655/CIRCexplorer2/known_circ.txt",
        find_circ="/home/choukaihsuan/GSE55872_results/circRNA/SRR444655/find_circ/splice_sites.bed",
        gt="/home/choukaihsuan/GSE113230_results/benchmark/rnaser/ground_truth.tsv",
    ),
    "NCI-H23": dict(
        ciriquant_gtf="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/circRNA/SRR17235469/SRR17235469.gtf",
        dcc="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/circRNA/SRR17235469/DCC/CircCoordinates",
        ce2="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/circRNA/SRR17235469/CIRCexplorer2/known_circ.txt",
        find_circ="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/circRNA/SRR17235469/find_circ/splice_sites.bed",
        gt="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/rnaser/ground_truth.tsv",
    ),
    "SW480": dict(
        ciriquant_gtf="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/circRNA/SRR17235468/SRR17235468.gtf",
        dcc="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/circRNA/SRR17235468/DCC/CircCoordinates",
        ce2="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/circRNA/SRR17235468/CIRCexplorer2/known_circ.txt",
        find_circ="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/circRNA/SRR17235468/find_circ/splice_sites.bed",
        gt="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/rnaser/ground_truth.tsv",
    ),
}


def build_index(coords):
    by_chr = {}
    for (c, s, e) in coords:
        by_chr.setdefault(c, []).append((s, e))
    for c in by_chr:
        by_chr[c].sort()
    return by_chr


def fast_match(chrom, start, end, idx, slop):
    lst = idx.get(chrom)
    if not lst:
        return False
    starts = [s for s, e in lst]
    lo = bisect.bisect_left(starts, start - slop)
    hi = bisect.bisect_right(starts, start + slop)
    for i in range(lo, hi):
        s2, e2 = lst[i]
        if max(abs(s2 - start), abs(e2 - end)) <= slop:
            return True
    return False


def load_ground_truth(path):
    gt = pd.read_csv(path, sep="\t")
    gt_pos, gt_neg = set(), set()
    for _, row in gt.iterrows():
        m = re.match(r'^(.+):(\d+)\|(\d+)$', row["circ_id"])
        if not m:
            continue
        chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        if row["is_true"] == 1:
            gt_pos.add((chrom, s, e))
        elif row["is_true"] == 0:
            gt_neg.add((chrom, s, e))
    return gt_pos, gt_neg


for cell, paths in DATASETS.items():
    print(f"\n{'='*70}\n{cell}\n{'='*70}")

    gt_pos, gt_neg = load_ground_truth(paths["gt"])
    gt_pos_idx, gt_neg_idx = build_index(gt_pos), build_index(gt_neg)
    print(f"ground truth: positive={len(gt_pos)} negative={len(gt_neg)}")

    ciri_map = cf.parse_ciriquant(paths["ciriquant_gtf"], MIN_BSJ, max_junction_ratio=1.0,
                                   qc_bsj_threshold=5)
    is_tp = {}
    for coord in ciri_map:
        chrom, start, end = coord
        pos_hit = fast_match(chrom, start, end, gt_pos_idx, SLOP_GT)
        neg_hit = fast_match(chrom, start, end, gt_neg_idx, SLOP_GT)
        if pos_hit and not neg_hit:
            is_tp[coord] = True
        elif neg_hit and not pos_hit:
            is_tp[coord] = False
        elif pos_hit and neg_hit:
            is_tp[coord] = True  # matches both -- same tie-break as stage0_partner_lift.py
        # else: not evaluable, excluded (no key added)

    n_evaluable = len(is_tp)
    overall_tp_rate = sum(is_tp.values()) / n_evaluable if n_evaluable else float("nan")
    print(f"CIRIquant calls (post pseudo-circ QC): {len(ciri_map)}; "
          f"ground-truth-evaluable: {n_evaluable}; overall TP rate: {overall_tp_rate:.4f}")

    # Partner coordinate sets, each tool's own slop recorded explicitly.
    partner_idx = {}
    for name, path, slop in [("DCC", paths["dcc"], PARTNER_SLOP["DCC"]),
                              ("CIRCexplorer2", paths["ce2"], PARTNER_SLOP["CIRCexplorer2"]),
                              ("find_circ", paths["find_circ"], PARTNER_SLOP["find_circ"])]:
        try:
            if name == "DCC":
                d = pd.read_csv(path, sep="\t")
                coords = {(r["Chr"], int(r["Start"]), int(r["End"])) for _, r in d.iterrows()}
            elif name == "CIRCexplorer2":
                coords = set(cf.parse_circexplorer2(path, MIN_BSJ).keys())
            else:
                coords = set(cf.parse_find_circ(path, MIN_BSJ).keys())
            partner_idx[name] = (build_index(coords), slop)
            print(f"  partner {name}: {len(coords)} calls, slop={slop}bp")
        except FileNotFoundError:
            print(f"  partner {name}: FILE NOT FOUND at {path} -- reporting gap, not substituting")

    # Stratify each evaluable CIRIquant call by how many partners corroborate it.
    support_counts = {}
    for coord in is_tp:
        chrom, start, end = coord
        n_support = 0
        for name, (idx, slop) in partner_idx.items():
            if fast_match(chrom, start, end, idx, slop):
                n_support += 1
        support_counts[coord] = n_support

    print(f"\n{'support_n':>10} {'n_loci':>8} {'n_TP':>8} {'TP_rate':>10} {'lift':>10}")
    for n in range(len(partner_idx) + 1):
        coords_at_n = [c for c, n_s in support_counts.items() if n_s == n]
        n_loci = len(coords_at_n)
        n_tp = sum(1 for c in coords_at_n if is_tp[c])
        tp_rate = n_tp / n_loci if n_loci else float("nan")
        lift = tp_rate / overall_tp_rate if n_loci and overall_tp_rate else float("nan")
        print(f"{n:>10} {n_loci:>8} {n_tp:>8} {tp_rate:>10.4f} {lift:>10.4f}")

print("\nDone. No conclusions written -- numbers only.")
