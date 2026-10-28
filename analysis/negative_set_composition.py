#!/usr/bin/env python3
"""
Negative-set BSJ support-read composition, full set (not just circDEX-only
false positives) across Hs68/NCI-H23/SW480. Read-only: only reads existing
ground_truth.tsv and DCC CircRNACount files. No Snakemake, no recomputation
of any detection/consensus step.

Reuses the exact DCC lookup convention already established and verified in
round4_b1_mechanism.py: CircRNACount's 4th column ("Chimeric.out.junction",
index 3) is the BSJ support-read count; a ground-truth locus is matched to
a DCC row within slop=10 bp (circDEX's own consensus slop), ties broken by
nearest coordinate distance -- same tie-break consensus_filter.py itself
uses when more than one candidate falls within slop.
"""
import bisect
import json
import sys
import pandas as pd
from scipy.stats import chi2_contingency

SLOP = 10
MIN_BSJ = 2

DATASETS = {
    "Hs68": dict(
        gt="/home/choukaihsuan/GSE113230_results/benchmark/rnaser/ground_truth.tsv",
        dcc="/home/choukaihsuan/GSE55872_results/circRNA/SRR444655/DCC/CircRNACount",
    ),
    "NCI-H23": dict(
        gt="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/rnaser/ground_truth.tsv",
        dcc="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_ncih23/circRNA/SRR17235469/DCC/CircRNACount",
    ),
    "SW480": dict(
        gt="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/rnaser/ground_truth.tsv",
        dcc="/home3/choukaihsuan/circRNA_agent/results/benchmark_multi/prjna789110_sw480/circRNA/SRR17235468/DCC/CircRNACount",
    ),
}

# The thesis's own §3.3 "circDEX-only negative calls" (b-set) exactly-2 figures,
# from audit/round4_audit.md (B1-3, verified against the corrected b-sets).
PAPER_B_SUBSET = {
    "Hs68":    (16, 20),    # not currently cited in thesis text
    "NCI-H23": (178, 217),
    "SW480":   (291, 363),
}


def load_dcc_counts(path):
    df = pd.read_csv(path, sep="\t")
    count_col = df.columns[3]
    by_chr = {}
    for _, r in df.iterrows():
        by_chr.setdefault(r["Chr"], []).append((int(r["Start"]), int(r["End"]), int(r[count_col])))
    for c in by_chr:
        by_chr[c].sort()
    return by_chr, count_col


def dcc_count_near(chrom, start, end, dcc_idx, slop):
    lst = dcc_idx.get(chrom)
    if not lst:
        return None
    starts = [s for s, e, c in lst]
    lo = bisect.bisect_left(starts, start - slop)
    hi = bisect.bisect_right(starts, start + slop)
    best, best_dist = None, None
    for i in range(lo, hi):
        s2, e2, c2 = lst[i]
        if max(abs(s2 - start), abs(e2 - end)) <= slop:
            dist = abs(s2 - start) + abs(e2 - end)
            if best is None or dist < best_dist:
                best, best_dist = c2, dist
    return best


def bucket(cnt):
    if cnt is None:
        return "no_dcc_match"
    if cnt == 2:
        return "eq2"
    if cnt == 3:
        return "eq3"
    if 4 <= cnt <= 10:
        return "4_10"
    if cnt > 10:
        return "gt10"
    return "lt_min_bsj"  # cnt < 2 -- shouldn't pass consensus, flagged if it appears


BUCKET_ORDER = ["eq2", "eq3", "4_10", "gt10", "no_dcc_match", "lt_min_bsj"]

results = {}
for cell, paths in DATASETS.items():
    truth = pd.read_csv(paths["gt"], sep="\t")
    dcc_idx, count_col = load_dcc_counts(paths["dcc"])
    gt_neg = truth[truth["is_true"] == 0]

    counts = {b: 0 for b in BUCKET_ORDER}
    raw_vals = []
    for _, row in gt_neg.iterrows():
        cid = row["circ_id"]
        chrom, rest = cid.split(":", 1)
        s, e = rest.split("|")
        s, e = int(s), int(e)
        cnt = dcc_count_near(chrom, s, e, dcc_idx, SLOP)
        counts[bucket(cnt)] += 1
        raw_vals.append(cnt)

    n_total = len(gt_neg)
    results[cell] = {"n_total": n_total, "counts": counts, "count_col": count_col}

    print(f"=== {cell} ===  N_neg_total = {n_total}  (DCC count field: '{count_col}')")
    for b in BUCKET_ORDER:
        n = counts[b]
        pct = n / n_total * 100 if n_total else float("nan")
        print(f"  {b:14s}: {n:6d}  ({pct:5.1f}%)")
    print()

# ── Cross-cell-line homogeneity test ────────────────────────────────────────
# Chi-square test of independence, not KS: KS is built for two continuous
# (or ordinal) sample CDFs; here we have THREE groups and a small number of
# discrete, unordered-for-the-test-purpose categories (bucket membership),
# which is exactly the chi-square-of-independence use case (row=cell line,
# col=bucket, H0: bucket proportions are independent of cell line). KS's
# precondition (two samples, continuous/ordinal comparison) does not hold
# for a 3-group categorical comparison. chi-square's own precondition
# (expected cell counts not too small, conventionally >=5 in at least 80% of
# cells) is checked below for each table rather than assumed.
#
# Two separate tables, answering two different questions:
#   (A) full composition INCLUDING no_dcc_match -- "is the complete negative
#       set's breakdown (has-a-DCC-call-at-all, and how strong) the same
#       across cell lines?"
#   (B) composition CONDITIONAL on having any DCC match at all -- "given a
#       negative locus that DCC did call, is the distribution of how
#       strongly it called it the same across cell lines?" This is the
#       question actually comparable to the paper's b-subset figures, since
#       the b-subset is by construction DCC-matched (circDEX's consensus
#       requires CIRIquant+DCC agreement).
cell_lines = list(DATASETS.keys())

print("=== Table A: full composition, including no_dcc_match ===")
buckets_A = ["eq2", "eq3", "4_10", "gt10", "no_dcc_match"]
table_A = [[results[cell]["counts"][b] for b in buckets_A] for cell in cell_lines]
print("          " + "  ".join(f"{b:>12s}" for b in buckets_A))
for cell, row in zip(cell_lines, table_A):
    print(f"{cell:10s}" + "  ".join(f"{v:12d}" for v in row))
chi2_A, p_A, dof_A, exp_A = chi2_contingency(table_A)
n_lt5_A = (exp_A < 5).sum()
print(f"chi2={chi2_A:.3f}  dof={dof_A}  p={p_A:.4g}  "
      f"({n_lt5_A}/{exp_A.size} expected cells < 5)")

print("\n=== Table B: composition conditional on DCC match (excludes no_dcc_match) ===")
buckets_B = ["eq2", "eq3", "4_10", "gt10"]
table_B = [[results[cell]["counts"][b] for b in buckets_B] for cell in cell_lines]
print("          " + "  ".join(f"{b:>8s}" for b in buckets_B))
for cell, row in zip(cell_lines, table_B):
    print(f"{cell:10s}" + "  ".join(f"{v:8d}" for v in row))
chi2_B, p_B, dof_B, exp_B = chi2_contingency(table_B)
n_lt5_B = (exp_B < 5).sum()
print(f"chi2={chi2_B:.3f}  dof={dof_B}  p={p_B:.4g}  "
      f"({n_lt5_B}/{exp_B.size} expected cells < 5, "
      f"min expected={exp_B.min():.2f})")

# ── Compare against the paper's circDEX-only "b"-subset figures ────────────
# Two comparisons, population-matched vs. not:
#   (1) full-set denominator (ALL negatives, mostly no-DCC-match) vs.
#       b-subset (ALWAYS DCC-matched, by construction) -- population
#       mismatch flagged explicitly, not silently presented as comparable.
#   (2) DCC-matched-subset denominator (excludes no_dcc_match) vs. b-subset
#       -- the population-matched comparison, since both sides are
#       conditioned on "DCC called this locus at all."
print("\n=== Comparison 1: FULL negative set (denominator includes no-DCC-match) vs. paper's b-subset ===")
print("(population mismatch: b-subset is by construction DCC-matched; full set mostly is not)")
for cell in cell_lines:
    full_eq2 = results[cell]["counts"]["eq2"]
    full_total = results[cell]["n_total"]
    full_pct = full_eq2 / full_total * 100
    b_eq2, b_total = PAPER_B_SUBSET[cell]
    b_pct = b_eq2 / b_total * 100
    print(f"{cell:10s} full-set exactly-2: {full_eq2}/{full_total} ({full_pct:.1f}%)   "
          f"b-subset (paper): {b_eq2}/{b_total} ({b_pct:.1f}%)   "
          f"delta: {full_pct - b_pct:+.1f} pp")

print("\n=== Comparison 2: DCC-MATCHED subset of full negative set vs. paper's b-subset (population-matched) ===")
for cell in cell_lines:
    c = results[cell]["counts"]
    matched_total = c["eq2"] + c["eq3"] + c["4_10"] + c["gt10"]
    matched_eq2 = c["eq2"]
    matched_pct = matched_eq2 / matched_total * 100 if matched_total else float("nan")
    b_eq2, b_total = PAPER_B_SUBSET[cell]
    b_pct = b_eq2 / b_total * 100
    print(f"{cell:10s} DCC-matched exactly-2: {matched_eq2}/{matched_total} ({matched_pct:.1f}%)   "
          f"b-subset (paper): {b_eq2}/{b_total} ({b_pct:.1f}%)   "
          f"delta: {matched_pct - b_pct:+.1f} pp")

print("\nDone. No conclusions written -- numbers only, per instruction.")
