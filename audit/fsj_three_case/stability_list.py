"""Effect-size stability and significance-retention for an explicit candidate list (edgeR only).

Input TSV columns: label, dataset, circ_id. For each row, the circRNA must be significant on the
full data (pvalue < 0.05 and |log2FC| >= 1) in that dataset to be tested; otherwise it is reported
as NOT_SIG_FULL and not tested.

Per run (one per dropped sample s), with full log2FC f and run log2FC l_s:
  unstable_effect = any(sign(l_s) != sign(f)) or any(|l_s| < 0.5 |f|)      (absent runs excluded)
  sig_dependent   = any(s: circRNA absent OR not (p<0.05 and |l_s|>=1))  (reported alongside)

usage: python stability_list.py <candidates.tsv> <out.tsv> [loo_root=~/loo_v2]
LOO root: ~/loo_v2 for every dataset (edgeR runs read from <root>/<GSE>/<drop>/edger.tsv).
"""
import os
import sys
import pandas as pd

CAND, OUT = sys.argv[1], sys.argv[2]
# LOO root: ~/loo_v2 (scripts/loo_stability.R; paired design read from patient_id, validated
# against the published tables for all 15 datasets). Override with a 3rd argument.
LOO_ROOT = sys.argv[3] if len(sys.argv) > 3 else os.path.expanduser("~/loo_v2")


def root(g):
    return LOO_ROOT


def load(g, s):
    fname = "edger.tsv"
    return pd.read_csv(os.path.join(root(g), g, s, fname), sep="\t").set_index("circ_id")


cand = pd.read_csv(CAND, sep="\t")
rows = []
for r in cand.itertuples():
    g, cid = r.dataset, r.circ_id
    full = load(g, "NONE")
    if cid not in full.index:
        rows.append(dict(label=r.label, dataset=g, circ_id=cid, status="NOT_TESTED_FULL")); continue
    f = full.loc[cid]
    fl, fp = float(f.log2FC), float(f.pvalue)
    if not (fp < 0.05 and abs(fl) >= 1):
        rows.append(dict(label=r.label, dataset=g, circ_id=cid, status="NOT_SIG_FULL",
                         log2FC_full=round(fl, 3), p_full=fp)); continue
    samples = sorted(x for x in os.listdir(os.path.join(root(g), g))
                       if x != "NONE" and os.path.isdir(os.path.join(root(g), g, x)))
    flips = drop50 = absent = sigloss = 0
    ratios = []
    lfcs = []
    for s in samples:
        t = load(g, s)
        if cid not in t.index:
            absent += 1; sigloss += 1; continue
        l = t.loc[cid]
        lv, pv = float(l.log2FC), float(l.pvalue)
        lfcs.append(round(lv, 2))
        if lv * fl < 0: flips += 1
        if abs(lv) < 0.5 * abs(fl): drop50 += 1
        ratios.append(abs(lv) / abs(fl))
        if not (pv < 0.05 and abs(lv) >= 1): sigloss += 1
    rows.append(dict(label=r.label, dataset=g, circ_id=cid, status="TESTED",
                     log2FC_full=round(fl, 3), p_full=fp,
                     min_abs_ratio=round(min(ratios), 3) if ratios else None,
                     n_sign_flip=flips, n_drop_gt50pct=drop50, n_absent=absent,
                     unstable_effect=bool(flips > 0 or drop50 > 0),
                     sig_dependent=bool(sigloss > 0),
                     n_runs=len(samples), log2FC_by_drop=";".join(map(str, lfcs))))

out = pd.DataFrame(rows)
out.to_csv(OUT, sep="\t", index=False)
pd.set_option("display.width", 250)
print(out.drop(columns=["log2FC_by_drop"], errors="ignore").to_string(index=False))
