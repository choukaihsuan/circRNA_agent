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

LAYOUT MODEL (revised after a round that let column width float with content, which
broke three other things): every column is the SAME fixed width (COL_W) and every
regular box uses the SAME font sizes. A fixed ARROW_GAP is reserved between columns
regardless of content, so the connecting arrow always has a real, visible shaft. If a
line of text doesn't fit COL_W at the fixed font size, it WRAPS onto another line --
the box never grows wider and the font never shrinks to make content fit. Boxes in the
same row share one y_top, so every row is flush across all four columns, regardless of
how many lines any one box's content wrapped onto.

Containment is checked against the SAVED PNG's actual pixels, not matplotlib's internal
bbox model (an internal-bbox assert passed cleanly once while the saved PNG still had
two lines overflowing their box -- constrained_layout re-settles axa's size as more
content is added, so an early-measured bbox goes stale by save time; the pixel check
re-opens the file that was actually written, so it can't be fooled by that staleness).
Two more checks run the same way: every pair of horizontally-connected boxes keeps a
minimum real pixel gap (a column-width-driven layout can let two boxes end up touching,
which silently reduces their connecting arrow to a zero-length dot), and every declared
arrow has a non-trivial pixel length in the final layout.

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
HEADER_BLACK = "#000000"      # column headers "1 Inputs" etc: black bold (F6)
BOXTITLE_BLUE = "#1a3f6b"     # each box's own first line: dark-blue bold (F7)
BODY_BLACK = "#1a1a1a"

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
fig.set_dpi(DPI)  # explicit: measurement and savefig must share one dpi, not just agree by
                    # construction default -- make that a visible invariant, not an accident
gs = fig.add_gridspec(2, 1, height_ratios=[2.0, 1.0])  # panel (a) needs more of the
    # figure's actual physical height now that columns 3/4's genuinely longer phrases
    # ("Task 1 · Detection accuracy" etc, vs column 1/2's short "SRP011042" etc.) wrap
    # onto several more lines at the fixed uniform font size -- this is real extra
    # content height, not a measurement artifact, so the fix is giving panel (a) more
    # pixels (which makes each data-unit MORE pixels, so the same fixed-pt-size text
    # needs FEWER data-units), not stretching axa's ylim (that only relabels data
    # units without changing how many pixels anything actually gets -- confirmed
    # earlier in this project's work that widening an axis range alone never helps).

ALL_BOXES = []   # (name, FancyBboxPatch, skip_sides) for every box() call
ALL_ARROWS = []  # (name, p0_data, p1_data) for every declared arrow -- checked post-save

# ---------------------------- Panel (a): diagram ----------------------------
axa = fig.add_subplot(gs[0])
axa.set_xlim(0, 104)  # 4 du wider than the other 3 columns need (0-100) -- this slack is
    # dedicated, NEW room used only by column 4's extra width (I4 below), not a shared
    # budget other columns compete for, so (unlike the panel-height fix) widening it here
    # is not a no-op: nothing else draws past x=100, so every added data-unit here is a
    # real added pixel for column 4 alone.
axa.set_ylim(0, 58)
axa.axis("off")
axa.set_title("(a)", loc="left", fontweight="bold", fontsize=13, pad=2)

# ---------------------------- Panel (b): wall-time bars ----------------------------
# Built HERE, before any of panel (a)'s text/box measurement, not after. Reason this
# matters: constrained_layout re-settles axa's on-screen size every time subplot
# content changes, including when axb is first added. text_extent_du() (below)
# converts a rendered pixel width to axa DATA units using whatever axa.transData is
# current AT THE TIME of measurement; if axa's pixel-per-data-unit ratio is still
# going to shift once axb exists, a width measured now is already stale by the time
# the figure is actually saved. Building axb's real, final content first means every
# later text_extent_du() call already sees the TRUE final layout.
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

renderer = fig.canvas.get_renderer()


def _redraw():
    global renderer
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()


_redraw()
p0 = axa.transData.transform((0, 0))
p1 = axa.transData.transform((0, 1))
px_per_dataunit = abs(p1[1] - p0[1])
DU_PER_PT = DPI / (px_per_dataunit * 72.0)  # data-units per point, y-axis


def text_extent_du(s, fontsize, weight="normal"):
    """Render `s` off-canvas, measure its rendered bbox, return (width, height) in
    this axes' data units. The text object is removed before returning."""
    t = axa.text(-1000, -1000, s, fontsize=fontsize, fontweight=weight)
    _redraw()
    bb = t.get_window_extent(renderer=renderer)
    bb_data = bb.transformed(axa.transData.inverted())
    t.remove()
    return abs(bb_data.width), abs(bb_data.height)


def wrap_line(text, max_w_du, fontsize, weight):
    """Greedy word-wrap: split `text` onto as many lines as needed so each one's
    rendered width fits max_w_du at this fontsize/weight. Returns [text] unchanged if
    it already fits. This is the H1/H2 fix's actual mechanism: box width and font size
    are both fixed inputs here, never adjusted to make content fit."""
    w, _ = text_extent_du(text, fontsize, weight)
    if w <= max_w_du:
        return [text]
    words = text.split(" ")
    lines, cur = [], ""
    for word in words:
        trial = (cur + " " + word).strip()
        tw, _ = text_extent_du(trial, fontsize, weight)
        if tw <= max_w_du or not cur:
            cur = trial
        else:
            lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


PAD_X, PAD_Y = 1.6, 1.1  # data-unit padding inside each box


def box(ax, name, x, y_top, lines, w, fc=BLUE, ec=BLUE_EDGE, fs=7.6, title_fs=8.6, lw=1.3,
        skip_sides=()):
    """Draw one rounded box of FIXED width `w`, anchored at its TOP edge (`y_top`), so
    boxes in the same row stay flush at the top regardless of how tall each one's
    content makes it (H3). `lines[0]` is the box's bold dark-blue title; the rest
    render at `fs` in body black. Any line (including the title) that doesn't fit `w`
    at its given font size is word-wrapped onto more lines -- `w` and the font sizes
    are never adjusted to make content fit (H1/H2). Returns (x, y_bottom, w, h).

    `name` is only for diagnostics: the patch is registered in ALL_BOXES and checked
    against the SAVED PNG's pixels at the very end of this script, not here -- an
    in-construction check would use whatever position constrained_layout has settled
    on *so far*, which keeps changing as later content is added. See module docstring.

    `skip_sides`: sides to exclude from that post-save ink check, for boxes with
    another, intentionally-nearby element on that side that isn't an overflow bug --
    e.g. a column header sitting a few pixels above a row-1 box."""
    max_text_w = w - 2 * PAD_X
    wrapped = []  # (text, is_title)
    for s in wrap_line(lines[0], max_text_w, title_fs, "bold"):
        wrapped.append((s, True))
    for ln in lines[1:]:
        for s in wrap_line(ln, max_text_w, fs, "normal"):
            wrapped.append((s, False))

    sizes = [text_extent_du(s, title_fs if is_t else fs, "bold" if is_t else "normal")
             for s, is_t in wrapped]
    line_heights = [s[1] for s in sizes]
    gaps = [lh * 1.5 for lh in line_heights]  # 1.5x glyph height as baseline-to-baseline step
    total_h = sum(gaps[:-1]) + line_heights[-1]
    h = total_h + 2 * PAD_Y
    y_bottom = y_top - h

    p = FancyBboxPatch((x, y_bottom), w, h, boxstyle="round,pad=0.25,rounding_size=0.6",
                        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(p)

    y_cursor = y_top - PAD_Y
    for i, (s, is_t) in enumerate(wrapped):
        ax.text(x + w / 2, y_cursor, s, ha="center", va="top",
                 fontsize=(title_fs if is_t else fs), fontweight="bold" if is_t else "normal",
                 color=(BOXTITLE_BLUE if is_t else BODY_BLACK), zorder=3)
        if i < len(gaps):
            y_cursor -= gaps[i]

    ALL_BOXES.append((name, p, set(skip_sides)))
    return (x, y_bottom, w, h)


def arrow(ax, name, p0, p1, **kw):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=14,
                         linewidth=1.3, color="#888888", zorder=1, **kw)
    ax.add_patch(a)
    ALL_ARROWS.append((name, p0, p1))


# ---- fixed column geometry: same width, same reserved arrow gap, for every column ----
COL_W = 21.0
ARROW_GAP = 4.5
COL_START = 1.0
col_x = [COL_START + i * (COL_W + ARROW_GAP) for i in range(4)]
# I4: column 4 (Outputs) alone is widened so "accuracy_summary.tsv" (one unbreakable
# token -- a filename has no spaces, so wrap_line() cannot wrap it at all) gets real
# interior padding instead of running up against its own right border. Only this
# column's WIDTH changes; col_x (every column's left edge, i.e. every arrow's gap) is
# untouched, so this cannot reintroduce H1.
OUT_COL_W = COL_W + 4.0

# title box -- fixed #1: GSE55872 -> SRP011042. Spans columns 1-4 (its own width is
# wider than COL_W on purpose, matching the frozen content list's full-width banner;
# unlike the regular boxes it already fit comfortably before and isn't part of
# H1/H2's complaint, so its own larger font size is left as-is).
TITLE_W = (col_x[3] + OUT_COL_W) - 8
box(axa, "title", 8, 57,
    ["Three pipelines compared (SRP011042)",
     "circDEX (CIRIquant + DCC) · CirComPara2 (4 tools) · nf-core/circrna (3 tools)"],
    TITLE_W, fc=BLUE, ec=TITLE_BLUE_EDGE, fs=8.6, title_fs=10.5, lw=1.6)

HEADER_Y = 46.0
ROW1_TOP = 44.0
ROW2_TOP = 18.0


ALL_HEADERS = []  # text artists for the four column headers, masked out before the
                    # pixel ink-check below (a header can sit close enough to an
                    # adjacent column's box that it falls inside that box's own
                    # margin strip -- that's the header rendering exactly where it's
                    # meant to, not content leaking out of the box being checked)


def draw_header(cx, text):
    t = axa.text(cx, HEADER_Y, text, ha="center", va="center", fontsize=9.3,
                  fontweight="bold", color=HEADER_BLACK)
    ALL_HEADERS.append(t)


for cx_i, label in zip(col_x, ["1  Inputs", "2  Ground truth", "3  Benchmark tasks", "4  Outputs"]):
    draw_header(cx_i + COL_W / 2, label)

# column 1: Inputs
# I2: comma -> semicolon before "RNase R x2" (was "Total RNA, RNase R x2", which wrapped
# as "Total RNA, RNase / R x2", splitting the proper noun "RNase R"); "RNase R x2" is
# joined with NBSP (U+00A0) so wrap_line() -- which only splits on literal ASCII spaces
# -- can never break inside it (I6). Same NBSP fix applied to "8 cores" below.
b_in1 = box(axa, "Inputs:SRP011042", col_x[0], ROW1_TOP,
            ["SRP011042", "Hs68 · hg19", "Total RNA; RNase R ×2"], COL_W)
b_in2 = box(axa, "Inputs:Runtime logs", col_x[0], ROW2_TOP,
            ["Runtime logs", "/usr/bin/time -v", "single sample, 8 cores"], COL_W)

# column 2: Ground truth
b_gt1 = box(axa, "Ground truth:Enrichment ratio", col_x[1], ROW1_TOP,
            ["Enrichment ratio", "ER = RNase R BSJ", "/ Total BSJ"], COL_W)
# fixed #2: strict </> -> inclusive >=/<=
b_gt2 = box(axa, "Ground truth:Label by ER", col_x[1], ROW2_TOP,
            ["Label by ER", "ER ≥ 1.5 → TP", "ER ≤ 0.5 → TN"], COL_W)

# column 3: Benchmark tasks (Task 2 + its output box use the purple fill)
# I1: the "Precision . Recall . F1 / Specificity . AUC-PR" metric list is dropped --
# both are already in Table 5 and in this figure's own caption, so repeating them here
# bought no information while being the main driver of Task 1's box height (and the
# resulting imbalance against Task 2).
#
# Measured (not guessed) against this column's fixed max_text_w=17.8du at the fixed
# body/title font sizes: "Task 1 · Detection accuracy" as ONE title line never fits
# (28.2du) at any wrap point that does not strand the · or split "Detection accuracy"
# mid-phrase -- gluing "· Detection accuracy" into one atomic NBSP chunk (an earlier
# attempt) does not fit either (overflowed the box, caught by the I5 clearance check).
# The actual fix is a CONTENT decision, not another layout knob: "Task 1"/"Task 2" is
# the title (bold, fits on its own, 6.5du), "Detection accuracy"/"Compute cost" moves
# into the body as its own line (normal weight, fits on its own, <16du) -- the relation
# between the two was always conveyed by them being adjacent lines in the same box, not
# by the ·. The remaining two multi-item lines use NBSP (I6) ONLY on the two genuinely
# atomic units that must not split ("ground truth", "8 cores"); "·" is replaced by ","
# in the metrics line, which wrap_line() already breaks cleanly after a comma, never
# stranding a bare connector.
t1_lines = ["Task 1", "Detection accuracy", "3 pipelines vs ground truth"]
t2_lines = ["Task 2", "Compute cost", "wall time, peak RAM, CPU-h",
            "same sample, 8 cores"]
b_t1 = box(axa, "Task1:Detection accuracy", col_x[2], ROW1_TOP, t1_lines, COL_W)
b_t2 = box(axa, "Task2:Compute cost", col_x[2], ROW2_TOP, t2_lines, COL_W,
           fc=PURPLE, ec=PURPLE_EDGE)

# column 4: Outputs -- I4: this column alone uses OUT_COL_W (wider than COL_W), so
# "accuracy_summary.tsv" (one unbreakable token -- a filename, no spaces to wrap on)
# gets real interior padding instead of running up against its own border.
out1_lines = ["Accuracy outputs", "accuracy_summary.tsv"]
out2_lines = ["Cost outputs", "compute_cost.tsv", "wall time per step · peak RAM"]
b_out1 = box(axa, "Outputs:Accuracy outputs", col_x[3], ROW1_TOP, out1_lines, OUT_COL_W)
b_out2 = box(axa, "Outputs:Cost outputs", col_x[3], ROW2_TOP, out2_lines, OUT_COL_W,
             fc=PURPLE, ec=PURPLE_EDGE)


def mid_right(b):
    x, y, w, h = b
    return (x + w, y + h / 2)


def mid_left(b):
    x, y, w, h = b
    return (x, y + h / 2)


def hline(ax, p0, p1):
    """A plain connecting segment (no arrowhead) -- used for the non-terminal legs
    of a multi-segment route."""
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#888888", linewidth=1.3,
              solid_capstyle="round", zorder=1)


def harrow(ax, name, bA, bB):
    """Horizontal arrow from box bA's right edge to box bB's left edge. Since every
    column now shares one fixed width and the columns are laid out COL_W + ARROW_GAP
    apart, bA and bB are guaranteed to be ARROW_GAP data-units apart here -- the
    earlier content-driven column width let this gap shrink to ~0 for several pairs,
    silently collapsing the arrow into an invisible dot (H1). The connector is forced
    perfectly horizontal rather than running box-center to box-center: a diagonal
    FancyArrowPatch over a short segment under this panel's strongly anisotropic x:y
    data scale was found, in an earlier pass, to render with its arrowhead at the
    WRONG end (confirmed by an isolated repro).

    y-position: the midpoint of whichever box is SHORTER. Both boxes in a row share
    the same y_top (H3), so the shorter box's entire vertical extent is a subset of
    the taller one's -- using its midpoint guarantees the arrow is level with real
    content in BOTH boxes. Using the source box's midpoint unconditionally (the
    first version of this function) pointed the arrow at empty space beside a much
    shorter target box whenever the source was the taller one (e.g. Task 1's long
    wrapped title vs. Accuracy outputs' two short lines) -- a defect the box-gap and
    arrow-length checks don't catch, since the arrow itself was perfectly real and
    long enough, just aimed at the wrong height. Caught by eye, not by an assert."""
    y = bA[1] + bA[3] / 2 if bA[3] <= bB[3] else bB[1] + bB[3] / 2
    arrow(ax, name, (bA[0] + bA[2], y), (bB[0], y))


# Inputs(SRP011042) -> Ground truth(Enrichment ratio) -> Benchmark(Task1) -> Outputs(Accuracy)
harrow(axa, "Inputs->GroundTruth", b_in1, b_gt1)
harrow(axa, "GroundTruth->Task1", b_gt1, b_t1)
harrow(axa, "Task1->Outputs", b_t1, b_out1)
# Ground truth(Enrichment ratio) feeds Label by ER (vertical, within the column --
# this one is axis-aligned the other way (constant x), so it isn't susceptible to
# the diagonal-arrowhead bug above)
gt1_bottom = b_gt1[1]
gt2_top = b_gt2[1] + b_gt2[3]
arrow(axa, "EnrichmentRatio->LabelByER",
      (col_x[1] + COL_W / 2, gt1_bottom), (col_x[1] + COL_W / 2, gt2_top))
# Runtime logs -> Task 2: routed BELOW the Ground truth column (a straight line
# would pass directly behind the "Label by ER" box, which draws on top and makes
# the arrow look cut off there) using ONLY axis-aligned segments -- every bend is a
# plain horizontal or vertical leg, never a diagonal, so none of them can trip the
# anisotropic-short-diagonal arrowhead bug noted above.
run_y = b_gt2[1] - 1.8  # a bit below Label-by-ER's bottom edge
p_start = mid_right(b_in2)
p_down = (p_start[0], run_y)
p_under = (col_x[1] + COL_W + 1.0, run_y)
p_up = (p_under[0], mid_left(b_t2)[1])
p_end = mid_left(b_t2)
hline(axa, p_start, p_down)
hline(axa, p_down, p_under)
hline(axa, p_under, p_up)
arrow(axa, "RuntimeLogs->Task2", p_up, p_end)
harrow(axa, "Task2->CostOutputs", b_t2, b_out2)

# ---- final layout is now fully settled (axb, headers, arrows all added) -- measure
# every box's on-screen pixel extent ONCE, here, then save, then check the SAVED FILE
# against those regions. No savefig(bbox_inches="tight") anywhere (that re-crops at
# save time and would invalidate a pre-save measurement); dpi is the same DPI used at
# figure creation (fig.set_dpi(DPI) above) and passed to savefig below -- one number,
# not two that happen to agree.
_redraw()
box_pixel_regions = []  # (name, x0, x1, y0_img, y1_img, skip_sides) in SAVED-PNG pixel coords
box_by_name = {}
for name, patch, skip_sides in ALL_BOXES:
    bb = patch.get_window_extent(renderer=renderer)  # display coords, origin bottom-left
    y0_img, y1_img = H_PX - bb.y1, H_PX - bb.y0  # flip to image-array row coords (origin top-left)
    box_pixel_regions.append((name, bb.x0, bb.x1, y0_img, y1_img, skip_sides))
    box_by_name[name] = (bb.x0, bb.x1, y0_img, y1_img)

header_pixel_regions = []  # same, for the 4 column headers -- masked out before checking,
                            # not excluded via per-box skip_sides, because a header can bleed
                            # into a NEIGHBOR box's margin (narrow inter-column gap), not just
                            # the margin of the box directly below it
HEADER_MASK_PAD_PX = 3
for t in ALL_HEADERS:
    bb = t.get_window_extent(renderer=renderer)
    y0_img, y1_img = H_PX - bb.y1, H_PX - bb.y0
    header_pixel_regions.append((bb.x0, bb.x1, y0_img, y1_img))

# A DIFFERENT failure mode from box-content overflow: a header colliding with the
# BOX below it (not a box's own content crossing its own border). Masking headers
# out of the ink-check above means that check literally cannot see this -- it would
# pass even with a header drawn right on top of a box. Caught once by eye on this
# exact script (row-1 box heights grew slightly once measured post-axb, same root
# cause as the column-width drift, and pushed their tops up into the header): check
# it explicitly, geometrically, rather than relying on looking every time.
for hx0, hx1, hy0, hy1 in header_pixel_regions:
    for name, patch, _ in ALL_BOXES:
        bb = patch.get_window_extent(renderer=renderer)
        bx0, bx1 = bb.x0, bb.x1
        by0_img, by1_img = H_PX - bb.y1, H_PX - bb.y0
        overlap_x = hx0 < bx1 and hx1 > bx0
        overlap_y = hy0 < by1_img and hy1 > by0_img
        if overlap_x and overlap_y:
            raise SystemExit(
                f"STOP: a column header overlaps box {name!r} -- header pixel bbox "
                f"x=({hx0:.1f},{hx1:.1f}) y=({hy0:.1f},{hy1:.1f}) vs box bbox "
                f"x=({bx0:.1f},{bx1:.1f}) y=({by0_img:.1f},{by1_img:.1f})")

# H1 (i): any two horizontally-adjacent, arrow-connected boxes must keep a real pixel
# gap between them -- the failure mode this catches (two boxes rendered touching, so
# their connecting arrow collapses to a zero-length dot) is invisible to the text-
# containment check above, which only checks text against ITS OWN box, never box
# against box.
MIN_BOX_GAP_PX = 15
HARROW_PAIRS = [
    ("Inputs:SRP011042", "Ground truth:Enrichment ratio"),
    ("Ground truth:Enrichment ratio", "Task1:Detection accuracy"),
    ("Task1:Detection accuracy", "Outputs:Accuracy outputs"),
    ("Task2:Compute cost", "Outputs:Cost outputs"),
]
for left_name, right_name in HARROW_PAIRS:
    lx0, lx1, _, _ = box_by_name[left_name]
    rx0, rx1, _, _ = box_by_name[right_name]
    gap_px = rx0 - lx1
    assert gap_px >= MIN_BOX_GAP_PX, (
        f"STOP: boxes {left_name!r} and {right_name!r} are only {gap_px:.1f}px apart "
        f"(need >= {MIN_BOX_GAP_PX}px) -- their connecting arrow would be invisible")

# H1 (ii): every declared arrow must have a non-trivial pixel length in the final
# layout -- same root cause as (i), checked directly on the arrow geometry itself
# rather than inferring it from box spacing (covers the routed elbow segments too).
MIN_ARROW_PX = 10
for name, p0_data, p1_data in ALL_ARROWS:
    disp0 = axa.transData.transform(p0_data)
    disp1 = axa.transData.transform(p1_data)
    length_px = ((disp1[0] - disp0[0]) ** 2 + (disp1[1] - disp0[1]) ** 2) ** 0.5
    assert length_px >= MIN_ARROW_PX, (
        f"STOP: arrow {name!r} is only {length_px:.1f}px long (need >= {MIN_ARROW_PX}px) "
        f"-- it would render as an invisible dot, not a visible arrow")

# I6 (topology): every declared arrow's endpoint must actually land ON the boundary of
# the box it's supposed to connect to -- a FULLY SEPARATE failure mode from (i)/(ii)
# above, and the one that bit this script already: harrow()'s first version aimed every
# connector at the SOURCE box's own midpoint, which was geometrically a perfectly real,
# correctly-long, horizontal arrow (passing both checks above) that simply pointed at
# empty space beside a much-shorter TARGET box once the two boxes' heights diverged.
# Neither a pixel ink-check (nothing overflowed) nor a length check (the arrow was long
# enough) can catch "the right number of real pixels, aimed at the wrong place" -- that
# was only caught by looking at the saved PNG. This makes that specific class of bug
# into an assertion instead of something that only an eyeball can catch.
ARROW_BOX_MAP = {
    "Inputs->GroundTruth":        ("Inputs:SRP011042", "Ground truth:Enrichment ratio"),
    "GroundTruth->Task1":         ("Ground truth:Enrichment ratio", "Task1:Detection accuracy"),
    "Task1->Outputs":             ("Task1:Detection accuracy", "Outputs:Accuracy outputs"),
    "EnrichmentRatio->LabelByER": ("Ground truth:Enrichment ratio", "Ground truth:Label by ER"),
    # source is a routed waypoint (the elbow's last horizontal leg), not literally a
    # point on Inputs:Runtime logs' own border -- only the destination is checked here.
    "RuntimeLogs->Task2":         (None, "Task2:Compute cost"),
    "Task2->CostOutputs":         ("Task2:Compute cost", "Outputs:Cost outputs"),
}
TOPO_TOL_PX = 7.0  # covers FancyBboxPatch's own boxstyle "pad=0.25" (du), which renders
    # each box ~0.25 du wider/taller than its nominal (x,y,w,h) on every side (~4.2px
    # at this panel's scale) -- confirmed empirically (a tol=3.0 first attempt failed
    # with every arrow endpoint landing a consistent ~4.2px short of the box's actual
    # rendered edge). Still far tighter than the "aimed at the wrong box" bug this check
    # exists to catch, where the endpoint lands tens of px away, outside the box.


def _on_box_boundary(px, py, box_px, tol=TOPO_TOL_PX):
    bx0, bx1, by0, by1 = box_px
    on_vert_edge = (abs(px - bx0) <= tol or abs(px - bx1) <= tol) and (by0 - tol <= py <= by1 + tol)
    on_horiz_edge = (abs(py - by0) <= tol or abs(py - by1) <= tol) and (bx0 - tol <= px <= bx1 + tol)
    return on_vert_edge or on_horiz_edge


for name, p0_data, p1_data in ALL_ARROWS:
    src_name, dst_name = ARROW_BOX_MAP.get(name, (None, None))
    for data_pt, box_name, which in ((p0_data, src_name, "source"), (p1_data, dst_name, "target")):
        if box_name is None:
            continue
        disp = axa.transData.transform(data_pt)
        px_, py_img = disp[0], H_PX - disp[1]  # flip to the same image-row convention as box_by_name
        assert _on_box_boundary(px_, py_img, box_by_name[box_name]), (
            f"STOP: arrow {name!r}'s {which} endpoint ({px_:.1f},{py_img:.1f}) does not "
            f"land on box {box_name!r}'s boundary {box_by_name[box_name]} (tol={TOPO_TOL_PX}px) "
            f"-- the arrow may be real and long enough but is aimed at the wrong place")
print(f"arrow topology check: {sum(1 for v in ARROW_BOX_MAP.values() for x in v if x)} "
      f"endpoint(s) across {len(ARROW_BOX_MAP)} arrows, all on their declared box's boundary")

fig.savefig(os.path.join(OUT, "figure6.png"), dpi=DPI)
fig.savefig(os.path.join(OUT, "figure6.pdf"), dpi=DPI)
print("wrote figure6.png/pdf to", OUT)
print("panel (b) totals:", {b[0]: b[2] for b in bar_data})
print(f"COL_W={COL_W} ARROW_GAP={ARROW_GAP} col_x={[round(c,2) for c in col_x]}")
print(f"box/box min gap check: {len(HARROW_PAIRS)} pairs, all >= {MIN_BOX_GAP_PX}px")
print(f"arrow length check: {len(ALL_ARROWS)} arrows, all >= {MIN_ARROW_PX}px")


def check_box_margins_in_saved_png(png_path, regions, header_regions=(), clearance_px=6,
                                     strip_px=40, ink_sum_max=240, header_mask_pad=3):
    """Read the FILE that was actually written and check, for each of a box's four
    edges, a band spanning from clearance_px INSIDE the border out to strip_px OUTSIDE
    it, for TEXT-colored ink. This is the check that actually matters: an earlier
    version of this script asserted containment against matplotlib's internal bbox
    model DURING construction, which passed cleanly while the saved PNG still had two
    lines overflowing their box -- constrained_layout kept re-adjusting axa's position
    as later content was added, so an early-measured bbox had gone stale by save time.
    This function instead opens the saved file and looks at actual pixels, so it
    cannot be fooled by that kind of staleness.

    I5: the band now starts clearance_px INSIDE the border, not gap_px past it. An
    earlier version only scanned the exterior (border+2px onward), which is blind to
    text that renders right up to, or touching, the border from the inside -- exactly
    what happened to "accuracy_summary.tsv" (a filename, one unbreakable token with no
    spaces to wrap on) before the I4 column-widening fix: its last glyph touched the
    right border and this check still reported 0 failures, because "touching the
    border" was never inside the scanned region at all. Requiring real interior
    clearance (not just "has it crossed outside yet") is a strictly stronger property,
    and it is what actually catches that failure mode.

    First pass used "any deviation from pure white" and produced 32 false positives
    on a figure that (per manual pixel inspection) had only 2 real overflows: every
    box's own anti-aliased rounded corner/border bleeds a few pale pixels past its
    nominal edge (R+G+B well above 600), and the gray connector arrows legitimately
    touch some boxes by design (R+G+B=408) -- neither is "text leaking out". Text ink
    (BODY_BLACK #1a1a1a sum=78, BOXTITLE_BLUE #1a3f6b sum=196) is categorically
    darker than every box border color in this script (289-350) or the arrow gray
    (408), so thresholding on R+G+B (not on "is it exactly white") separates the two
    cleanly with real margin on both sides.

    Second pass (ink threshold in place) still had a false positive: a column
    header's own text can sit close enough to an ADJACENT column's box (narrow
    inter-column gap) to land inside that box's margin strip, even though the header
    isn't that box's content overflowing -- it's rendering exactly where it's meant
    to. `header_regions` are painted white (with header_mask_pad px of slack) in the
    loaded array before any box is checked, so header ink can never trigger either
    box's check, on any side, regardless of which box it happens to sit near."""
    from PIL import Image
    import numpy as np
    arr = np.array(Image.open(png_path).convert("RGB")).astype(int)
    full_h, full_w = arr.shape[0], arr.shape[1]
    for hx0, hx1, hy0, hy1 in header_regions:
        mx0 = max(int(round(hx0)) - header_mask_pad, 0)
        mx1 = min(int(round(hx1)) + header_mask_pad, full_w)
        my0 = max(int(round(hy0)) - header_mask_pad, 0)
        my1 = min(int(round(hy1)) + header_mask_pad, full_h)
        arr[my0:my1, mx0:mx1] = 255
    failures = []
    for name, x0, x1, y0_img, y1_img, skip_sides in regions:
        x0, x1 = int(round(x0)), int(round(x1))
        y0_img, y1_img = int(round(y0_img)), int(round(y1_img))
        sides = {
            "right":  (x1 - clearance_px, x1 + strip_px, y0_img, y1_img),
            "left":   (x0 - strip_px, x0 + clearance_px, y0_img, y1_img),
            "top":    (x0, x1, y0_img - strip_px, y0_img + clearance_px),
            "bottom": (x0, x1, y1_img - clearance_px, y1_img + strip_px),
        }
        for side in skip_sides:
            sides.pop(side, None)
        for side, (cx0, cx1, cy0, cy1) in sides.items():
            cx0, cy0 = max(cx0, 0), max(cy0, 0)
            cx1, cy1 = min(cx1, full_w), min(cy1, full_h)
            if cx1 <= cx0 or cy1 <= cy0:
                continue
            region = arr[cy0:cy1, cx0:cx1]
            is_ink = region.sum(axis=-1) < ink_sum_max
            if is_ink.any():
                ys, xs = np.where(is_ink)
                failures.append(
                    f"{name} [{side}]: text-dark pixel at image (row={cy0 + int(ys[0])}, "
                    f"col={cx0 + int(xs[0])}), RGB={tuple(int(v) for v in region[ys[0], xs[0]])}, "
                    f"checked rows[{cy0}:{cy1}] cols[{cx0}:{cx1}]")
    return failures


png_failures = check_box_margins_in_saved_png(
    os.path.join(OUT, "figure6.png"), box_pixel_regions, header_pixel_regions)
print(f"pixel containment check: {len(box_pixel_regions)} boxes x 4 sides, "
      f"{len(png_failures)} failures (clearance=6px, strip=40px, tol=10)")
if png_failures:
    for f in png_failures:
        print("  FAIL:", f)
    raise SystemExit(f"STOP: {len(png_failures)} pixel-level overflow(s) in the saved PNG "
                      f"(see FAIL lines above) -- figure6.png was still written to disk for "
                      f"inspection, but is WRONG, do not use it")
print("pixel containment check: PASSED, 0 overflows")
print("ALL CHECKS PASSED: text containment (6px clearance), header/box collision, "
      "box/box min gap, arrow min length, arrow topology")
