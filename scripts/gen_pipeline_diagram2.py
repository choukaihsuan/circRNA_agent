"""
Generate CircDEX Analysis Pipeline diagram v2 (pipeline_diagram2.png/svg).
Changes from v1:
  - Title capitalised to "CircDEX Analysis Pipeline"
  - Download box expanded: shows fallback chain + fasterq-dump + pigz compression
  - Download box taller to accommodate three lines of sub-text
Run: conda run -n ciriquant python scripts/gen_pipeline_diagram2.py [output_dir]
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patches as mpatches
import sys, os

FIG_W, FIG_H = 14.4, 7.6
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, FIG_W)
ax.set_ylim(0, FIG_H)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── Palette ───────────────────────────────────────────────
C_IN    = '#1e3a5f'
C_QC    = '#0b7864'
C_DET   = '#4a1c87'
C_FILT  = '#7a1414'
C_ANNOT = '#786600'
C_RPT   = '#1a5c20'
C_ARROW = '#6b7280'

# ── Y positions ───────────────────────────────────────────
Y_TITLE = 7.44
Y_CIRI  = 6.48
Y_STAR  = 5.55
Y_MAIN  = 5.88   # GEO / Download / QC
Y2      = 3.78
Y3      = 2.35
Y4      = 0.97
Y_LEG   = 0.36

# ── Arrow / line helpers ──────────────────────────────────
ASTYLE = mpatches.ArrowStyle('->', head_width=0.16, head_length=0.14)

def arrow(x1, y1, x2, y2):
    p = FancyArrowPatch((x1, y1), (x2, y2),
                        arrowstyle=ASTYLE, color=C_ARROW,
                        linewidth=1.8, zorder=4,
                        connectionstyle='arc3,rad=0')
    ax.add_patch(p)

def hline(x1, y, x2):
    ax.plot([x1, x2], [y, y], color=C_ARROW, lw=1.8, zorder=2, solid_capstyle='butt')

def vline(x, y1, y2):
    ax.plot([x, x], [y1, y2], color=C_ARROW, lw=1.8, zorder=2, solid_capstyle='butt')

# ── Box helper ────────────────────────────────────────────
def rbox(cx, cy, w, h, color, title, sub='', ts=9.5, ss=7.8, lh=0.19, z=3):
    r = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                       boxstyle='round,pad=0.07',
                       linewidth=0, facecolor=color, zorder=z)
    ax.add_patch(r)
    if sub:
        ax.text(cx, cy + lh, title, ha='center', va='center',
                color='white', fontsize=ts, fontweight='bold', zorder=z+1)
        ax.text(cx, cy - lh * 0.85, sub, ha='center', va='center',
                color='white', fontsize=ss, zorder=z+1, linespacing=1.32)
    else:
        ax.text(cx, cy, title, ha='center', va='center',
                color='white', fontsize=ts, fontweight='bold', zorder=z+1)

# ── Detection group background ────────────────────────────
ax.add_patch(FancyBboxPatch((5.18, 5.00), 5.40, 2.05,
                            boxstyle='round,pad=0.14',
                            linewidth=1.2, edgecolor='#c4b5e8',
                            facecolor='#f3f0f8', zorder=1))

# ── Title ─────────────────────────────────────────────────
ax.text(FIG_W/2, Y_TITLE, 'CircDEX Analysis Pipeline',
        ha='center', va='center', fontsize=17, fontweight='bold', color='#111827')

# ── Boxes ─────────────────────────────────────────────────
# Row 1 left
rbox(0.92,  Y_MAIN, 1.40, 0.88, C_IN, 'GEO / SRA', 'raw data', ts=9, ss=8.0)

# Download — expanded: 3-line sub (fallback chain + fasterq-dump + pigz)
rbox(2.82,  Y_MAIN, 1.88, 1.04, C_QC, 'Download',
     '① aria2c S3  ② ascp\n③ prefetch (HTTPS)\nfasterq-dump → pigz',
     ts=9.5, ss=7.3, lh=0.24)

rbox(4.64,  Y_MAIN, 1.56, 0.88, C_QC, 'QC / Trim',
     'FastQC · fastp ·\nMultiQC', ss=7.8)

# Row 1 right (Detection)
rbox(7.38, Y_CIRI, 2.32, 0.73, C_DET, 'CIRIquant',
     'HISAT2 + BWA-MEM → GTF', ss=8.2, lh=0.17)
rbox(7.38, Y_STAR, 2.32, 0.75, C_DET, 'STAR × 3',
     'paired · mate1 · mate2', ss=8.2, lh=0.17)
rbox(10.40, Y_STAR, 1.92, 0.75, C_DET, 'DCC 0.5.0',
     'CircCoordinates', ss=8.2, lh=0.17)

# Row 2
rbox(1.28, Y2, 2.22, 0.87, C_FILT, 'Consensus Filter',
     'min_tools=2 · slop=10bp\n· pseudo-circ QC', ts=9, ss=7.6, lh=0.19)
rbox(3.60, Y2, 1.74, 0.87, C_QC,   'Merge Counts', 'BSJ · FSJ matrix', lh=0.18)
rbox(5.52, Y2 + 0.245, 1.60, 0.43, C_ANNOT, 'Assign Isoforms',    ts=9)
rbox(5.52, Y2 - 0.245, 1.60, 0.43, C_ANNOT, 'circBase Annotation', ts=9)
rbox(9.06, Y2, 4.58, 0.93, C_FILT, 'DE Analysis (×3)',
     'edgeR_ciriquant · DESeq2 · limma-voom', ts=10, ss=8.5, lh=0.21)

# Row 3
rbox(1.84, Y3, 3.22, 0.87, C_IN, 'Predict Interactions',
     'CircInteractome · ENCORI · union top-50', ss=7.8, lh=0.19)
rbox(5.52, Y3, 3.72, 0.87, C_IN, 'Isoform Switching',
     'IUI · Wilcoxon · within-gene BH FDR', ss=7.8, lh=0.19)
rbox(9.78, Y3, 3.88, 0.87, C_IN, 'Rank Biomarkers',
     '6D score: sig · FC · conf · circBase · miRNA · RBP', ss=7.5, lh=0.19)

# Row 4
rbox(7.08, Y4, 13.50, 0.93, C_RPT, 'report.html',
     'Plotly · Type I/II · Venn diagram · SVG exon diagram · Biomarker table · Isoform switching',
     ts=13, ss=9.0, lh=0.24)

# ── Arrows ────────────────────────────────────────────────
# ① GEO → Download → QC
arrow(1.62, Y_MAIN, 1.88, Y_MAIN)
DL_R = 2.82 + 1.88/2    # 3.76
QC_L = 4.64 - 1.56/2    # 3.86
arrow(DL_R, Y_MAIN, QC_L, Y_MAIN)

# ② QC → branch → CIRIquant + STAR
QC_R = 4.64 + 1.56/2    # 5.42
BPX  = 5.85
hline(QC_R, Y_MAIN, BPX)
vline(BPX, Y_STAR, Y_CIRI)
hline(BPX, Y_CIRI, 6.22)
arrow(6.22, Y_CIRI, 6.24, Y_CIRI)
hline(BPX, Y_STAR, 6.22)
arrow(6.22, Y_STAR, 6.24, Y_STAR)

# ③ STAR → DCC
arrow(8.54, Y_STAR, 9.44, Y_STAR)

# ④ CIRIquant + DCC → Consensus Filter
Y_COL = Y2 + 0.87/2 + 0.20
CIRI_X = 7.38
DCC_X  = 10.40 + 1.92/2
CONS_X = 1.28

vline(CIRI_X, Y_CIRI - 0.73/2, Y_COL)
vline(DCC_X,  Y_STAR - 0.75/2, Y_COL)
hline(CONS_X, Y_COL, DCC_X)
vline(CONS_X, Y2 + 0.87/2, Y_COL)
arrow(CONS_X, Y2 + 0.87/2, CONS_X, Y2 + 0.87/2 - 0.02)

# ⑤ Consensus → Merge
arrow(1.28 + 2.22/2, Y2, 2.73, Y2)

# ⑥ Merge → Assign Isoforms + circBase Annotation
MERGE_R = 3.60 + 1.74/2
FORK_X  = 4.73
ANNO_L  = 5.52 - 1.60/2

hline(MERGE_R, Y2, FORK_X)
vline(FORK_X, Y2 - 0.245, Y2 + 0.245)
hline(FORK_X, Y2 + 0.245, ANNO_L)
arrow(ANNO_L, Y2 + 0.245, ANNO_L + 0.02, Y2 + 0.245)
hline(FORK_X, Y2 - 0.245, ANNO_L)
arrow(ANNO_L, Y2 - 0.245, ANNO_L + 0.02, Y2 - 0.245)

# ⑦ Annotations → DE Analysis
COLL2_X = 6.52
ANNO_R  = 5.52 + 1.60/2
DE_L    = 9.06 - 4.58/2

hline(ANNO_R, Y2 + 0.245, COLL2_X)
hline(ANNO_R, Y2 - 0.245, COLL2_X)
vline(COLL2_X, Y2 - 0.245, Y2 + 0.245)
hline(COLL2_X, Y2, DE_L)
arrow(DE_L, Y2, DE_L + 0.02, Y2)

# ⑧ DE → three output nodes
DE_CX  = 9.06
DE_BOT = Y2 - 0.93/2
Y_FAN  = Y3 + 0.87/2 + 0.12

for tx in [1.84, 5.52, 9.78]:
    vline(DE_CX, DE_BOT, Y_FAN)
    hline(DE_CX, Y_FAN, tx)
    vline(tx, Y_FAN, Y3 + 0.87/2)
    arrow(tx, Y3 + 0.87/2, tx, Y3 + 0.87/2 - 0.02)

# ⑨ Three output nodes → report.html
RPT_TOP = Y4 + 0.93/2
for sx in [1.84, 5.52, 9.78]:
    vline(sx, Y3 - 0.87/2, RPT_TOP)
    arrow(sx, RPT_TOP, sx, RPT_TOP - 0.02)

# ── Legend ────────────────────────────────────────────────
items = [
    (C_IN,    'Input / Interaction'),
    (C_QC,    'QC / Count'),
    (C_DET,   'Detection'),
    (C_FILT,  'Filter / DE'),
    (C_ANNOT, 'Annotation'),
    (C_RPT,   'Report output'),
]
LX0, LBW, LBH, GAP = 0.55, 0.38, 0.30, 2.20
for i, (col, lbl) in enumerate(items):
    lx = LX0 + i * GAP
    ax.add_patch(FancyBboxPatch((lx, Y_LEG - LBH/2), LBW, LBH,
                                boxstyle='round,pad=0.04',
                                linewidth=0, facecolor=col, zorder=3))
    ax.text(lx + LBW + 0.13, Y_LEG, lbl,
            ha='left', va='center', fontsize=8.5, color='#374151')

# ── Save ──────────────────────────────────────────────────
out_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
plt.tight_layout(pad=0)
for ext in ('png', 'svg'):
    path = os.path.join(out_dir, f'pipeline_diagram2.{ext}')
    fig.savefig(path, dpi=130, bbox_inches='tight', facecolor='white',
                format=ext if ext == 'svg' else None)
print(f"Saved → {out_dir}/pipeline_diagram2.{{png,svg}}")
