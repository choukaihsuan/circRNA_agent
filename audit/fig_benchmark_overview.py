"""Figure 6: benchmark overview (panel (a) pipeline/task diagram + panel (b) wall-time bars).

This replaces a previously undocumented, unreproducible PNG (word/media/image10.png in
the thesis docx) that could not be traced to any generating script despite a thorough
search (local repo, server repo, Claude Projects folder, byte-size/pixel-dimension match
against every PNG reachable) -- see audit/ for the provenance trail. That image also had
two now-confirmed text errors, which this script corrects at their only two true sources
(panel (a) is where both strings live -- fixing them anywhere else, e.g. in a report HTML
or a different pptx, does not touch this figure):

  1. "GSE55872" -> "SRP011042" (title box + Inputs box): GSE55872 is Drosophila in GEO;
     the dataset actually used is SRP011042 (Jeck et al. 2013, Hs68 human fibroblasts).
  2. "ER > 1.5 -> TP" / "ER < 0.5 -> TN" -> "ER >= 1.5" / "ER <= 0.5": the ground-truth
     script (rnaser_ground_truth.py) uses inclusive boundaries; strict "</>" does not
     reproduce the published 3,221/1,664/291 positive/negative/excluded split.

Everything else -- box text, box order, box colors, arrow routing, the plain bold
Arabic 1/2/3/4 column numbers (not circled), panel (b)'s layout -- is kept identical to
the frozen content list the thesis author read off the published image. Panel (b)'s
numbers are read from compute_cost.tsv at run time, never hardcoded.

usage: fig_benchmark_overview.py <compute_cost.tsv> <out_dir>
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
})
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

TSV, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

# ---- target pixel size: matches the superseded image10.png exactly (2177x2074) ----
DPI = 300
W_PX, H_PX = 2177, 2074
FIG_W, FIG_H = W_PX / DPI, H_PX / DPI

BLUE = "#d9e8f8"
BLUE_EDGE = "#4472a8"
PURPLE = "#e6e0f2"
PURPLE_EDGE = "#6b5b95"
TITLE_BLUE_EDGE = "#2c5c99"

# ============================== Panel (a): data ==============================

def rd_row(rows, pipeline):
    for r in rows:
        if r["Pipeline"] == pipeline:
            return r
    raise SystemExit(f"STOP: no compute_cost.tsv row for pipeline {pipeline!r}")


def fnum(v):
    return float(v) if v not in (None, "") else 0.0


with open(TSV, newline="") as fh:
    cost_rows = list(csv.DictReader(fh, delimiter="\t"))

TOOLS = ["CIRIquant_min", "STAR_min", "DCC_min", "CIRCexplorer2_min", "find_circ_min"]
TOOL_LABELS = ["CIRIquant", "STAR", "DCC", "CIRCexplorer2", "find_circ"]
TOOL_COLORS = ["#1f3f66", "#4a79b5", "#b9d3ea", "#e2792e", "#5fa15a"]

BAR_ROWS = [
    ("circDEX", "Our pipeline"),
    ("CirComPara2", "CirComPara2_4tools"),
    ("nf-core/circrna", "nf-core/circrna"),
]
bar_data = []
for label, pid in BAR_ROWS:
    row = rd_row(cost_rows, pid)
    vals = [fnum(row[t]) for t in TOOLS]
    total = fnum(row["Total_wall_min"])
    bar_data.append((label, vals, total))

# ============================== figure / layout ==============================

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI, constrained_layout=True)
gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0])

# ---------------------------- Panel (a): diagram ----------------------------
axa = fig.add_subplot(gs[0])
axa.set_xlim(0, 100)
axa.set_ylim(0, 58)
axa.axis("off")
axa.set_title("(a)", loc="left", fontweight="bold", fontsize=13, pad=2)

# Force constrained_layout to settle axa's final on-screen position *before* we place
# any text, so we can convert "points" -> "data units" exactly instead of guessing a
# gap and hoping it doesn't overlap (that guess was wrong the first time around: text
# anchors a few data-units apart looked fine on paper but the actual glyphs, at
# 9-10pt, were taller than that gap once rendered).
fig.canvas.draw()
p0 = axa.transData.transform((0, 0))
p1 = axa.transData.transform((0, 1))
px_per_dataunit = abs(p1[1] - p0[1])
DU_PER_PT = DPI / (px_per_dataunit * 72.0)  # data-units per point, y-axis


def line_gap(fontsize, leading=1.35):
    """Vertical data-unit offset between successive text baselines at this fontsize."""
    return fontsize * leading * DU_PER_PT


def box(ax, x, y, w, h, lines, fc=BLUE, ec=BLUE_EDGE, fs=7.6, title_fs=None, lw=1.3):
    """Draw one rounded box, vertically centered text stack. If title_fs is given,
    `lines[0]` is rendered bold at title_fs and the rest at fs; otherwise all lines
    render at fs (no bold first line)."""
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.6",
                        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(p)
    cy = y + h / 2
    if title_fs:
        gaps = [line_gap(title_fs)] + [line_gap(fs)] * (len(lines) - 2)
    else:
        gaps = [line_gap(fs)] * (len(lines) - 1)
    total_h = sum(gaps)
    y_cursor = cy + total_h / 2
    for i, ln in enumerate(lines):
        is_title = title_fs and i == 0
        ax.text(x + w / 2, y_cursor, ln, ha="center", va="center",
                 fontsize=(title_fs if is_title else fs),
                 fontweight="bold" if is_title else "normal",
                 color="#1a1a1a", zorder=3)
        if i < len(gaps):
            y_cursor -= gaps[i]
    return (x, y, w, h)


def arrow(ax, p0, p1, **kw):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=14,
                         linewidth=1.3, color="#888888", zorder=1, **kw)
    ax.add_patch(a)


# title box -- fixed #1: GSE55872 -> SRP011042
box(axa, 8, 50, 84, 7.0,
    ["Three pipelines compared (SRP011042)",
     "circDEX (CIRIquant + DCC) · CirComPara2 (4 tools) · nf-core/circrna (3 tools)"],
    fc=BLUE, ec=TITLE_BLUE_EDGE, fs=8.6, title_fs=10.5, lw=1.6)

col_x = [1, 26.5, 52, 77.5]
col_w = 21
headers = ["1  Inputs", "2  Ground truth", "3  Benchmark tasks", "4  Outputs"]
for cx, hd in zip(col_x, headers):
    axa.text(cx + col_w / 2, 44.6, hd, ha="center", va="center", fontsize=9.3,
              fontweight="bold", color="#17315c")

# column 1: Inputs
b_in1 = box(axa, col_x[0], 34.0, col_w, 8.8,
            ["SRP011042", "Hs68 · hg19", "Total RNA, RNase R ×2"], fs=7.6)
b_in2 = box(axa, col_x[0], 22.0, col_w, 7.3,
            ["Runtime logs", "/usr/bin/time -v", "single sample, 8 cores"], fs=7.6)

# column 2: Ground truth
b_gt1 = box(axa, col_x[1], 34.0, col_w, 8.8,
            ["Enrichment ratio", "ER = RNase R BSJ", "/ Total BSJ"], fs=7.6)
# fixed #2: strict </> -> inclusive >=/<=
b_gt2 = box(axa, col_x[1], 22.0, col_w, 7.3,
            ["Label by ER", "ER ≥ 1.5 → TP", "ER ≤ 0.5 → TN"], fs=7.6)

# column 3: Benchmark tasks (Task 2 + its output box use the purple fill)
b_t1 = box(axa, col_x[2], 31.5, col_w, 11.5,
           ["Task 1 · Detection accuracy", "3 pipelines vs ground truth",
            "Precision · Recall · F1", "Specificity · AUC-PR"],
           fs=7.1, title_fs=7.6)
b_t2 = box(axa, col_x[2], 19.5, col_w, 8.8,
           ["Task 2 · Compute cost", "wall time · peak RAM · CPU-h",
            "same sample, 8 cores"], fc=PURPLE, ec=PURPLE_EDGE, fs=7.1, title_fs=7.6)

# column 4: Outputs
b_out1 = box(axa, col_x[3], 31.5, col_w, 11.5,
             ["Accuracy outputs", "accuracy_summary.tsv"], fs=7.3, title_fs=7.8)
b_out2 = box(axa, col_x[3], 19.5, col_w, 8.8,
             ["Cost outputs", "compute_cost.tsv", "wall time per step · peak RAM"],
             fc=PURPLE, ec=PURPLE_EDGE, fs=7.1, title_fs=7.6)


def mid_right(b):
    x, y, w, h = b
    return (x + w, y + h / 2)


def mid_left(b):
    x, y, w, h = b
    return (x, y + h / 2)


# Inputs(SRP011042) -> Ground truth(Enrichment ratio) -> Benchmark(Task1) -> Outputs(Accuracy)
arrow(axa, mid_right(b_in1), mid_left(b_gt1))
arrow(axa, mid_right(b_gt1), mid_left(b_t1))
arrow(axa, mid_right(b_t1), mid_left(b_out1))
# Ground truth(Enrichment ratio) feeds Label by ER (vertical, within the column)
arrow(axa, (col_x[1] + col_w / 2, 34.0), (col_x[1] + col_w / 2, 29.3))
# Runtime logs -> Task 2: routed BELOW the Ground truth column (a straight line here
# would pass directly behind/through the "Label by ER" box, which draws over it and
# makes the arrow look like it dead-ends there instead of reaching Task 2).
dip_x, dip_y = col_x[1] + col_w / 2, 17.0
start = mid_right(b_in2)
axa.plot([start[0], dip_x], [start[1], dip_y], color="#888888", linewidth=1.3,
          solid_capstyle="round", zorder=1)
arrow(axa, (dip_x, dip_y), mid_left(b_t2))
arrow(axa, mid_right(b_t2), mid_left(b_out2))

# ---------------------------- Panel (b): wall-time bars ----------------------------
axb = fig.add_subplot(gs[1])
axb.set_title("(b)", loc="left", fontweight="bold", fontsize=13, pad=2)

y_pos = list(range(len(bar_data)))[::-1]
for (label, vals, total), y in zip(bar_data, y_pos):
    left = 0.0
    for v, c in zip(vals, TOOL_COLORS):
        if v > 0:
            axb.barh(y, v, left=left, color=c, edgecolor="none", height=0.55)
        left += v
    axb.text(left + 15, y, f"{total:,.0f} min", va="center", ha="left",
              fontsize=10, color="#555555")

axb.set_yticks(y_pos)
axb.set_yticklabels([b[0] for b in bar_data], fontsize=11)
axb.set_xlim(0, 1400)
axb.set_xlabel("Wall time (min · single sample · 8 cores)", fontsize=10)
axb.set_ylim(-0.7, len(bar_data) - 0.3)
axb.spines[["top", "right", "left"]].set_visible(False)
axb.grid(axis="x", color="#e5e7eb", linewidth=0.8, zorder=0)
axb.set_axisbelow(True)
axb.tick_params(axis="y", length=0)

handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in TOOL_COLORS]
axb.legend(handles, TOOL_LABELS, loc="upper center", bbox_to_anchor=(0.5, -0.18),
           ncol=5, frameon=False, fontsize=9.5, handlelength=1.2, handleheight=1.2)

fig.savefig(os.path.join(OUT, "figure6.png"), dpi=DPI)
fig.savefig(os.path.join(OUT, "figure6.pdf"), dpi=DPI)
print("wrote figure6.png/pdf to", OUT)
print("panel (b) totals:", {b[0]: b[2] for b in bar_data})
print("DU_PER_PT (sanity check):", DU_PER_PT)
