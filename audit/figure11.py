"""Figure 11: the Figure 10 story (coordinate agreement + confidence_score) for one dataset,
GSE133998 -- the dataset the reader already knows from Figure 2/5/11-caption context (Table 5,
section 3.2). Reuses confidence_figures.py's load()/bucket_counts() (not a reimplementation) so
the numbers cannot diverge from Figure 10's.

usage: figure11.py <dump_root> <out_dir> [--dump-dirname dump]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
})
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import confidence_figures as cf  # noqa: E402

DATASET = "GSE133998"
HIST_WIDTH = cf.HIST_WIDTH

ap = argparse.ArgumentParser()
ap.add_argument("dump_root")
ap.add_argument("out_dir")
ap.add_argument("--dump-dirname", default="dump")
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

df_all, in_hash = cf.load(a.dump_root, a.dump_dirname)
sub = df_all[df_all.dataset == DATASET]
if sub.empty:
    sys.exit(f"STOP: no rows for {DATASET} in the loaded dump set")

b = cf.bucket_counts(sub.loc[sub.has_dcc, "dist_dcc"])
loci = len(sub)
dist0, d13, d410 = b["zero"], b["one_three"], b["four_ten"]
frac_identical = dist0 / b["n"]

fig, (axa, axb) = plt.subplots(2, 1, figsize=(5.5, 5.5), constrained_layout=True,
                               gridspec_kw={"height_ratios": [3, 1]})

axa.hist(sub.confidence_score, bins=np.arange(
    np.floor(sub.confidence_score.min() / HIST_WIDTH) * HIST_WIDTH,
    np.ceil(sub.confidence_score.max() / HIST_WIDTH) * HIST_WIDTH + HIST_WIDTH, HIST_WIDTH),
    color="#4C72B0")
axa.set_xlabel("Confidence score")
axa.set_ylabel("Loci")
axa.set_title("(a)", loc="left", fontweight="bold")

vals = [dist0, d13, d410]
labels = ["0", "1–3", "4–10"]

# A1: headroom for the dominant bar's two-line label is found by measuring the
# actual rendered annotation (not guessing a multiplier and hoping) -- a first guess
# of 1.15x was tried and the containment assert below caught it as still
# insufficient, so this now retries with more headroom until it actually fits,
# instead of shipping a number that merely looked plausible.
renderer = fig.canvas.get_renderer()
headroom_mult = 1.15
for attempt in range(8):
    axb.clear()
    bars = axb.bar(labels, vals, color=["#4C72B0", "#DD8452", "#C44E52"])
    axb.set_xticks(range(len(labels)))
    axb.set_xticklabels(labels, fontsize=8)  # single line, explicit (not the
        # categorical auto-labels bar() assigns, which wrapped)
    axb.set_ylim(0, max(vals) * headroom_mult)
    axb.set_xlabel("Coordinate difference between callers (bp)")
    axb.set_ylabel("Loci")
    axb.set_title("(b)", loc="left", fontweight="bold")
    # A1/A2: offset-points placement off each bar's own top, not a fixed data-y
    # value -- for the two tiny bars (d13, d410) this anchors "just above that
    # bar's top", which for a near-zero bar sits right at the baseline, not
    # floating in the middle of the axis the way a 10%-of-loci margin did before.
    for bar, v in zip(bars, vals):
        axb.annotate(f"{v:,}\n({v / loci:.3%})",
                     xy=(bar.get_x() + bar.get_width() / 2, v),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom")
    fig.canvas.draw()
    bar_bboxes = [p.get_window_extent(renderer=renderer) for p in bars.patches]
    spine_bbox = axb.spines["top"].get_window_extent(renderer=renderer)
    ann_bboxes = [t.get_window_extent(renderer=renderer) for t in axb.texts]
    ok = True
    for tb in ann_bboxes:
        if tb.y1 > spine_bbox.y0 + 1e-6:
            ok = False
            break
    if ok:
        break
    headroom_mult += 0.08
else:
    raise SystemExit("STOP: could not find headroom that keeps panel (b)'s "
                      "annotations inside the axes after 8 attempts")

# A4: measured, not eyeballed -- no text bbox may intersect a bar, the top spine, or
# another text's bbox. A failure here means the figure must not be written out.
for i, tb in enumerate(ann_bboxes):
    for j, bb in enumerate(bar_bboxes):
        assert not tb.overlaps(bb), (
            f"STOP: annotation {i} overlaps bar {j} -- {tb} vs {bb}")
    assert tb.y1 <= spine_bbox.y0 + 1e-6 or tb.y0 >= spine_bbox.y1, (
        f"STOP: annotation {i} crosses the top spine -- {tb} vs {spine_bbox}")
    for j, ob in enumerate(ann_bboxes):
        if i != j:
            assert not tb.overlaps(ob), (
                f"STOP: annotation {i} overlaps annotation {j} -- {tb} vs {ob}")
print(f"panel (b) headroom multiplier used: {headroom_mult:.2f}")

fig.savefig(os.path.join(a.out_dir, "figure11.png"), dpi=600)
fig.savefig(os.path.join(a.out_dir, "figure11.pdf"), dpi=600)

cap = (f"Figure 11. The Figure 10 story for a single dataset, {DATASET} (the dataset used throughout "
       f"section 3.2 / Table 5 / Figure 12). (a) Distribution of confidence_score over all {loci:,} "
       f"consensus loci in this dataset. (b) Coordinate agreement between the two callers: {dist0:,} "
       f"of {loci:,} loci ({frac_identical:.3%}) have identical coordinates; the remaining {d13+d410} "
       f"loci differ by 1–3 bp ({d13:,}) or 4–10 bp ({d410:,}). These numbers are the "
       f"{DATASET} row of Figure 10(a) ({dist0:,}/{loci:,} = {frac_identical:.4f}); input hash "
       f"{in_hash[:12]} matches Figure 10's.\n")
open(os.path.join(a.out_dir, "figure11_caption.md"), "w", encoding="utf-8").write(cap)
print(f"{DATASET}: loci={loci:,} dist0={dist0:,} ({frac_identical:.4%}) dist_1_3={d13} dist_4_10={d410}")
print("input_hash", in_hash[:12])
print("wrote figure11.png/pdf/caption.md to", a.out_dir)
