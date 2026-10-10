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

Box sizes for columns 3/4 are computed from the *actual rendered* text extent
(fig.canvas.draw() + get_window_extent()), not a fixed guessed width.

Containment is checked against the SAVED PNG's actual pixels, not matplotlib's internal
bbox model. An internal-bbox assert was tried first and it passed cleanly while the saved
PNG still had two lines overflowing their box -- the failure mode is specifically
"the internal model says it's fine, the rendered output isn't", which checking the
internal model can never catch, because constrained_layout keeps re-adjusting axa's
position/size as more content (axb, headers added after a box's own containment was
checked) gets added, so an early-measured bbox goes stale by save time. The pixel check
below runs once, at the very end, against the file that was actually written.

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
gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0])

ALL_BOXES = []  # (name, FancyBboxPatch) for every box() call, checked post-save at the end

# ---------------------------- Panel (a): diagram ----------------------------
axa = fig.add_subplot(gs[0])
axa.set_xlim(0, 100)
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
# the figure is actually saved -- which is exactly what happened on the previous
# pass: col3_w/col4_w were computed before axb existed, two box titles ended up
# wider than their (too-narrow) box in the saved PNG despite the internal numbers
# looking consistent at computation time. Building axb's real, final content first
# means every later text_extent_du() call already sees the TRUE final layout.
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


PAD_X, PAD_Y = 1.6, 1.1  # data-unit padding inside each box


def box(ax, name, x, y_bottom, lines, fc=BLUE, ec=BLUE_EDGE, fs=7.6, title_fs=9.0, lw=1.3,
        min_w=None, skip_sides=()):
    """Draw one rounded box sized to fit its own text exactly (measured, not guessed).
    `lines[0]` is the box's bold dark-blue title line; the rest render at `fs` in
    body black. Box width = widest rendered line + padding (or `min_w` if larger, so
    columns can still share a width when desired). Returns (x, y_bottom, w, h).

    `name` is only for diagnostics: the patch is registered in ALL_BOXES and checked
    against the SAVED PNG's pixels at the very end of this script (see
    check_box_margins_in_saved_png), not here -- an in-construction check here would
    use whatever position constrained_layout has settled on *so far*, which keeps
    changing as later content (more boxes, axb) is added, so it can pass while the
    final saved file still overflows. See the module docstring.

    `skip_sides`: sides to exclude from that post-save check, for boxes with another,
    intentionally-nearby element on that side that isn't an overflow bug -- e.g. the
    four top-row boxes each sit a few pixels below their own column's bold header
    ("1 Inputs" etc.), so a "top" check there trips on that header's own ink, not on
    this box's content leaking out."""
    sizes = [text_extent_du(lines[0], title_fs, "bold")]
    for ln in lines[1:]:
        sizes.append(text_extent_du(ln, fs, "normal"))
    w = max(s[0] for s in sizes) + 2 * PAD_X
    if min_w:
        w = max(w, min_w)
    line_heights = [s[1] for s in sizes]
    gaps = [h * 1.5 for h in line_heights]  # 1.5x glyph height as baseline-to-baseline step
    total_h = sum(gaps[:-1]) + line_heights[-1]
    h = total_h + 2 * PAD_Y

    p = FancyBboxPatch((x, y_bottom), w, h, boxstyle="round,pad=0.25,rounding_size=0.6",
                        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2)
    ax.add_patch(p)

    cy_top = y_bottom + h - PAD_Y  # baseline of the first (title) line
    y_cursor = cy_top
    text_objs = []
    for i, ln in enumerate(lines):
        is_title = i == 0
        t = ax.text(x + w / 2, y_cursor, ln, ha="center", va="top",
                     fontsize=(title_fs if is_title else fs),
                     fontweight="bold" if is_title else "normal",
                     color=(BOXTITLE_BLUE if is_title else BODY_BLACK), zorder=3)
        text_objs.append(t)
        if i < len(gaps):
            y_cursor -= gaps[i]

    ALL_BOXES.append((name, p, set(skip_sides)))  # checked against the saved PNG's pixels at the end
    return (x, y_bottom, w, h)


def arrow(ax, p0, p1, **kw):
    a = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=14,
                         linewidth=1.3, color="#888888", zorder=1, **kw)
    ax.add_patch(a)


# title box -- fixed #1: GSE55872 -> SRP011042
box(axa, "title", 8, 50,
    ["Three pipelines compared (SRP011042)",
     "circDEX (CIRIquant + DCC) · CirComPara2 (4 tools) · nf-core/circrna (3 tools)"],
    fc=BLUE, ec=TITLE_BLUE_EDGE, fs=8.6, title_fs=10.5, lw=1.6, min_w=84)

col_x = [1, 25, 49, None]  # col_x[3] is filled in once column 3's measured width is known
col_w_12 = 21  # columns 1-2 already fit at this fixed width; only 3-4 need measuring
HEADER_Y = 47.3


ALL_HEADERS = []  # text artists for the four column headers, masked out before the
                    # pixel ink-check below (a header can sit close enough to an
                    # adjacent column's box that it falls inside that box's own
                    # margin strip -- that's the header rendering exactly where it's
                    # meant to, not content leaking out of the box being checked)


def draw_header(cx, text):
    t = axa.text(cx, HEADER_Y, text, ha="center", va="center", fontsize=9.3,
                  fontweight="bold", color=HEADER_BLACK)
    ALL_HEADERS.append(t)


draw_header(col_x[0] + col_w_12 / 2, "1  Inputs")
draw_header(col_x[1] + col_w_12 / 2, "2  Ground truth")
# headers for columns 3 and 4 are drawn later, once their measured widths/positions
# are known -- drawing a guessed position first and painting over it with a second
# copy left a visible sliver of the first draw in an earlier version of this script.

# column 1: Inputs
b_in1 = box(axa, "Inputs:SRP011042", col_x[0], 34.0,
            ["SRP011042", "Hs68 · hg19", "Total RNA, RNase R ×2"],
            fs=7.6, title_fs=8.6, min_w=col_w_12)
b_in2 = box(axa, "Inputs:Runtime logs", col_x[0], 21.0,
            ["Runtime logs", "/usr/bin/time -v", "single sample, 8 cores"],
            fs=7.6, title_fs=8.6, min_w=col_w_12)

# column 2: Ground truth
b_gt1 = box(axa, "Ground truth:Enrichment ratio", col_x[1], 34.0,
            ["Enrichment ratio", "ER = RNase R BSJ", "/ Total BSJ"],
            fs=7.6, title_fs=8.6, min_w=col_w_12)
# fixed #2: strict </> -> inclusive >=/<=
b_gt2 = box(axa, "Ground truth:Label by ER", col_x[1], 21.0,
            ["Label by ER", "ER ≥ 1.5 → TP", "ER ≤ 0.5 → TN"],
            fs=7.6, title_fs=8.6, min_w=col_w_12)

# column 3: Benchmark tasks (Task 2 + its output box use the purple fill). Width is
# measured from content (previously a fixed guess caused the overflow in F1/F2) and
# shared between the two boxes in this column so their left/right edges align.
t1_lines = ["Task 1 · Detection accuracy", "3 pipelines vs ground truth",
            "Precision · Recall · F1", "Specificity · AUC-PR"]
t2_lines = ["Task 2 · Compute cost", "wall time · peak RAM · CPU-h",
            "same sample, 8 cores"]
T3_FS, T3_TITLE_FS = 5.3, 5.8
col3_w = max(
    max(text_extent_du(t1_lines[0], T3_TITLE_FS, "bold")[0],
        max(text_extent_du(l, T3_FS, "normal")[0] for l in t1_lines[1:])),
    max(text_extent_du(t2_lines[0], T3_TITLE_FS, "bold")[0],
        max(text_extent_du(l, T3_FS, "normal")[0] for l in t2_lines[1:])),
) + 2 * PAD_X
b_t1 = box(axa, "Task1:Detection accuracy", col_x[2], 31.0, t1_lines, fs=T3_FS,
           title_fs=T3_TITLE_FS, min_w=col3_w)
b_t2 = box(axa, "Task2:Compute cost", col_x[2], 18.5, t2_lines, fc=PURPLE, ec=PURPLE_EDGE,
           fs=T3_FS, title_fs=T3_TITLE_FS, min_w=col3_w)
draw_header(col_x[2] + col3_w / 2, "3  Benchmark tasks")

# column 4: Outputs -- same measure-first approach
out1_lines = ["Accuracy outputs", "accuracy_summary.tsv"]
out2_lines = ["Cost outputs", "compute_cost.tsv", "wall time per step · peak RAM"]
T4_FS, T4_TITLE_FS = 5.5, 5.9
col4_w = max(
    max(text_extent_du(out1_lines[0], T4_TITLE_FS, "bold")[0],
        max(text_extent_du(l, T4_FS, "normal")[0] for l in out1_lines[1:])),
    max(text_extent_du(out2_lines[0], T4_TITLE_FS, "bold")[0],
        max(text_extent_du(l, T4_FS, "normal")[0] for l in out2_lines[1:])),
) + 2 * PAD_X
col_x[3] = col_x[2] + col3_w + 0.5
b_out1 = box(axa, "Outputs:Accuracy outputs", col_x[3], 31.0, out1_lines, fs=T4_FS,
             title_fs=T4_TITLE_FS, min_w=col4_w)
b_out2 = box(axa, "Outputs:Cost outputs", col_x[3], 18.5, out2_lines, fc=PURPLE,
             ec=PURPLE_EDGE, fs=T4_FS, title_fs=T4_TITLE_FS, min_w=col4_w)

right_edge = col_x[3] + col4_w
xlim_right = axa.get_xlim()[1]
assert right_edge <= xlim_right, (
    f"STOP: column 4 content (right edge {right_edge:.2f}) overflows the panel's "
    f"xlim (0, {xlim_right}) -- widen the canvas or shrink columns 3/4 further, "
    f"don't just eyeball the render.")

draw_header(col_x[3] + col4_w / 2, "4  Outputs")


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


def harrow(ax, bA, bB):
    """Horizontal arrow from box bA's right edge to box bB's left edge, both ends at
    bA's own mid-height. Box heights are now sized to content, so two boxes in the
    same "row" rarely share an exact mid-y; a diagonal FancyArrowPatch between them
    (short segment + the panel's strongly anisotropic x:y data scale) was found to
    render with its arrowhead at the WRONG end (confirmed by an isolated repro) --
    forcing the connector perfectly horizontal (or, for the routed one below,
    perfectly horizontal/vertical per segment) avoids the diagonal entirely rather
    than trying to out-guess when the matplotlib bug does or doesn't trigger."""
    y = bA[1] + bA[3] / 2
    arrow(ax, (bA[0] + bA[2], y), (bB[0], y))


# Inputs(SRP011042) -> Ground truth(Enrichment ratio) -> Benchmark(Task1) -> Outputs(Accuracy)
harrow(axa, b_in1, b_gt1)
harrow(axa, b_gt1, b_t1)
harrow(axa, b_t1, b_out1)
# Ground truth(Enrichment ratio) feeds Label by ER (vertical, within the column --
# this one is axis-aligned the other way (constant x), so it isn't susceptible to
# the diagonal-arrowhead bug above)
gt1_bottom = b_gt1[1]
gt2_top = b_gt2[1] + b_gt2[3]
arrow(axa, (col_x[1] + col_w_12 / 2, gt1_bottom), (col_x[1] + col_w_12 / 2, gt2_top))
# Runtime logs -> Task 2: routed BELOW the Ground truth column (a straight line
# would pass directly behind the "Label by ER" box, which draws on top and makes
# the arrow look cut off there) using ONLY axis-aligned segments -- every bend is a
# plain horizontal or vertical leg, never a diagonal, so none of them can trip the
# anisotropic-short-diagonal arrowhead bug noted above.
run_y = b_gt2[1] - 1.8  # a bit below Label-by-ER's bottom edge
p_start = mid_right(b_in2)
p_down = (p_start[0], run_y)
p_under = (col_x[1] + col_w_12 + 1.0, run_y)
p_up = (p_under[0], mid_left(b_t2)[1])
p_end = mid_left(b_t2)
hline(axa, p_start, p_down)
hline(axa, p_down, p_under)
hline(axa, p_under, p_up)
arrow(axa, p_up, p_end)
harrow(axa, b_t2, b_out2)

# ---- final layout is now fully settled (axb, headers, arrows all added) -- measure
# every box's on-screen pixel extent ONCE, here, then save, then check the SAVED FILE
# against those regions. No savefig(bbox_inches="tight") anywhere (that re-crops at
# save time and would invalidate a pre-save measurement); dpi is the same DPI used at
# figure creation (fig.set_dpi(DPI) above) and passed to savefig below -- one number,
# not two that happen to agree.
_redraw()
box_pixel_regions = []  # (name, x0, x1, y0_img, y1_img, skip_sides) in SAVED-PNG pixel coords
for name, patch, skip_sides in ALL_BOXES:
    bb = patch.get_window_extent(renderer=renderer)  # display coords, origin bottom-left
    y0_img, y1_img = H_PX - bb.y1, H_PX - bb.y0  # flip to image-array row coords (origin top-left)
    box_pixel_regions.append((name, bb.x0, bb.x1, y0_img, y1_img, skip_sides))

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

fig.savefig(os.path.join(OUT, "figure6.png"), dpi=DPI)
fig.savefig(os.path.join(OUT, "figure6.pdf"), dpi=DPI)
print("wrote figure6.png/pdf to", OUT)
print("panel (b) totals:", {b[0]: b[2] for b in bar_data})
print("col3_w", col3_w, "col4_w", col4_w, "col_x[3]", col_x[3])


def check_box_margins_in_saved_png(png_path, regions, header_regions=(), gap_px=2,
                                     strip_px=8, ink_sum_max=240, header_mask_pad=3):
    """Read the FILE that was actually written and check a thin strip of pixels just
    outside each of a box's four edges (gap_px past the border, strip_px wide) for
    TEXT-colored ink. This is the check that actually matters: an earlier version of
    this script asserted containment against matplotlib's internal bbox model DURING
    construction, which passed cleanly while the saved PNG still had two lines
    overflowing their box -- constrained_layout kept re-adjusting axa's position as
    later content was added, so an early-measured bbox had gone stale by save time.
    This function instead opens the saved file and looks at actual pixels, so it
    cannot be fooled by that kind of staleness.

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
            "right":  (x1 + gap_px, x1 + gap_px + strip_px, y0_img, y1_img),
            "left":   (x0 - gap_px - strip_px, x0 - gap_px, y0_img, y1_img),
            "top":    (x0, x1, y0_img - gap_px - strip_px, y0_img - gap_px),
            "bottom": (x0, x1, y1_img + gap_px, y1_img + gap_px + strip_px),
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
      f"{len(png_failures)} failures (gap=2px, strip=8px, tol=10)")
if png_failures:
    for f in png_failures:
        print("  FAIL:", f)
    raise SystemExit(f"STOP: {len(png_failures)} pixel-level overflow(s) in the saved PNG "
                      f"(see FAIL lines above) -- figure6.png was still written to disk for "
                      f"inspection, but is WRONG, do not use it")
print("pixel containment check: PASSED, 0 overflows")
