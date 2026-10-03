"""Paper-style Figures 7, 8, 9 (15-dataset version only), styled after the v19 figures.
Outputs figs/paper_15ds/{fig7,fig8,fig9}.png (+pdf). Fig 9 is one 2737x3294 px PNG like the existing image13.png.
Colors: Fig7 #2E75B6/#E07B39; Fig8 #C0504D/#4E9A51; grid #E6E6E6; n.s. cell #EDEDED; Fig 9a bars and (b) number color
sampled from image13.png (#1F4E79 #2E75B6 #5B9BD5 #9DC3E6 #D9D9D9); heatmap scale = #C0504D - white - #4E9A51.
Run: python3 make_paper_figs_15ds.py   (matplotlib, numpy; inputs in ./data)
"""
import os, gzip, math, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from thesis_common import *
plt.rcParams.update({"font.family": "DejaVu Serif", "pdf.fonttype": 42})
summ, sig, tested, genes = load()
cb = {}
for line in gzip.open(os.path.join(DATA, "circbase_map.tsv.gz"), "rt"):
    a = line.rstrip("\n").split("\t")
    if len(a) == 2 and a[0] != "circ_id": cb[a[0]] = a[1]
GRID = "#E6E6E6"; SMALL_N = 20
NAME = {  # row/column labels in the v19 figure format "<type organ> (<accession>)"
 "GSE171011": "PTC thyroid", "GSE113230": "TNBC breast", "GSE136569": "PDAC pancreas", "SRP156355": "Breast IDC",
 "GSE97239": "Bladder", "PRJNA553289": "SCLC lung", "GSE133998": "Breast", "GSE121842": "CRC colorectal",
 "GSE221107": "Prostate", "GSE143797": "NPC nasopharyngeal", "GSE108735": "RCC kidney", "GSE248612": "Gastric",
 "GSE77509": "HCC liver", "GSE130078": "ESCC esophageal", "GSE148036": "LUAD lung"}
ORGAN_ORDER = {"bladder": 0, "breast": 1, "colorectal": 2, "esophagus": 3, "kidney": 4, "liver": 5, "lung": 6, "nasopharynx": 7, "ovary": 8, "pancreas": 9, "prostate": 10, "stomach": 11, "thyroid": 12}
def lab(d): return "%s  (%s)" % (NAME[d], d)
def nsig(d): return len(sig[d])
def up(d): return sum(1 for l, _ in sig[d].values() if l > 0)
OUT = os.path.join(HERE, "figs", "paper_15ds"); os.makedirs(OUT, exist_ok=True)
def save(fig, name): fig.savefig(os.path.join(OUT, name+".png"), dpi=200); fig.savefig(os.path.join(OUT, name+".pdf")); plt.close(fig)
def style(ax):
    ax.spines[["top", "right"]].set_visible(False)

# ---------------- Figure 7
def fig7():
    C1, C2 = "#2E75B6", "#E07B39"
    order = sorted(DS15, key=lambda d: -nsig(d)); n = len(order)
    fig, ax = plt.subplots(figsize=(10.5, 0.62*n+2.0)); ax.set_position([0.33, 0.14, 0.58, 0.84])
    ax.set_xticks(range(0, 101, 20)); ax.xaxis.grid(True, color=GRID, lw=1); ax.set_axisbelow(True)
    for i, d in enumerate(order):
        c = collections.Counter(t for _, t in sig[d].values()); a, b = c["Type_I"], c["Type_II"]; t = a+b; y = n-1-i; hp = "///" if t < SMALL_N else None
        ax.barh(y, 100*a/t, height=0.62, color=C1, hatch=hp, edgecolor="white", lw=0); ax.barh(y, 100*b/t, left=100*a/t, height=0.62, color=C2, hatch=hp, edgecolor="white", lw=0)
        ax.text(101.5, y, "n=%d%s" % (t, " †" if t < SMALL_N else ""), va="center", fontsize=12, fontweight="bold" if t < SMALL_N else "normal")
    ax.set_yticks(range(n)); ax.set_yticklabels([lab(d) for d in order[::-1]], fontsize=13); ax.set_xlim(0, 100); ax.set_ylim(-0.6, n-0.4)
    ax.set_xlabel("Share of significant DE circRNAs (%)", fontsize=14); ax.tick_params(axis="x", labelsize=12); style(ax)
    fig.legend([Patch(color=C1), Patch(color=C2), Patch(facecolor="#BBBBBB", edgecolor="white", hatch="///")],
               ["Type I  (circRNA-specific)", "Type II  (host-gene driven)", "† n < %d" % SMALL_N], loc="lower center", ncol=3, fontsize=12, frameon=False)
    save(fig, "fig7")

# ---------------- Figure 8
def fig8():
    CU, CD = "#C0504D", "#4E9A51"
    order = sorted(DS15, key=lambda d: (-up(d)/nsig(d), -nsig(d))); n = len(order)
    fig, ax = plt.subplots(figsize=(10.5, 0.62*n+2.0)); ax.set_position([0.33, 0.14, 0.52, 0.84])
    ax.set_xticks(range(0, 101, 20)); ax.xaxis.grid(True, color=GRID, lw=1); ax.set_axisbelow(True)
    for i, d in enumerate(order):
        u = up(d); dn = nsig(d)-u; t = nsig(d); y = n-1-i; pu = 100*u/t; note = []
        ax.barh(y, pu, height=0.62, color=CU, edgecolor="white", lw=0); ax.barh(y, 100-pu, left=pu, height=0.62, color=CD, edgecolor="white", lw=0)
        for val, x0, x1, col in ((u, 0, pu, CU), (dn, pu, 100, CD)):
            if val == 0:
                ax.text(2.4 if col == CU else 97.6, y, "0", ha="center", va="center", fontsize=11, color=col, fontweight="bold", bbox=dict(boxstyle="square,pad=0.2", fc="white", ec=col, lw=1.3))
            elif x1-x0 >= 7: ax.text((x0+x1)/2, y, str(val), ha="center", va="center", fontsize=11, color="white", fontweight="bold" if t < SMALL_N else "normal")
            else: note.append("%s%d" % ("↑" if col == CU else "↓", val))
        ax.text(101.5, y, "n=%d%s%s" % (t, " †" if t < SMALL_N else "", " (%s)" % note[0] if note else ""), va="center", fontsize=12, fontweight="bold" if t < SMALL_N else "normal")
    ax.set_yticks(range(n)); ax.set_yticklabels([lab(d) for d in order[::-1]], fontsize=13); ax.set_xlim(0, 100); ax.set_ylim(-0.6, n-0.4)
    ax.set_xlabel("Share of significant DE circRNAs (%)", fontsize=14); ax.tick_params(axis="x", labelsize=12); style(ax)
    fig.legend([Patch(color=CU), Patch(color=CD)], ["Significant — Up", "Significant — Down"], loc="lower center", ncol=2, fontsize=12, frameon=False)
    save(fig, "fig8")

# ---------------- Figure 9 (single PNG, 2737 x 3294 px)
def fig9():
    cnt = tiers(sig, DS15); mx = max(cnt.values())
    ge = {k: sum(1 for v in cnt.values() if v >= k) for k in range(1, mx+1)}
    fig = plt.figure(figsize=(13.685, 16.47), dpi=200)
    # (a)
    axa = fig.add_axes([0.14, 0.765, 0.30, 0.165]); levels = [5, 4, 3, 2, 1]
    vals = [sum(1 for v in cnt.values() if v == k) for k in levels]
    axa.bar(range(5), vals, color=["#1F4E79", "#2E75B6", "#5B9BD5", "#9DC3E6", "#D9D9D9"], width=0.8)
    axa.set_yscale("log"); axa.set_ylim(2, 3000); axa.set_xticks(range(5)); axa.set_xticklabels(levels, fontsize=14); style(axa)
    axa.tick_params(axis="y", labelsize=12); axa.set_xlabel("Datasets in which significant", fontsize=15); axa.set_ylabel("circRNAs (log)", fontsize=15)
    fig.text(0.115, 0.955, "(a)", fontsize=17, fontweight="bold")
    # (b)
    fig.text(0.52, 0.955, "(b)", fontsize=17, fontweight="bold")
    n5 = ge[5]
    for i, (num, txt) in enumerate([(ge[2], "significant in ≥ 2 datasets"), (ge[3], "recurrent in ≥ 3 datasets"), (ge[4], "recurrent in ≥ 4 datasets (%d of them in 5)" % n5)]):
        y = 0.918 - i*0.044
        fig.text(0.535, y, str(num), fontsize=32, fontweight="bold", color="#1F4E79", va="center")
        fig.text(0.635, y, txt, fontsize=15, va="center")
    # (c)
    rows = sorted([k for k, v in cnt.items() if v >= 4], key=lambda k: (genes.get(k, "~"), k))
    cols = sorted(DS15, key=lambda d: (ORGAN_ORDER[META[d][0]], NAME[d], d))
    M = np.full((len(rows), len(cols)), np.nan)
    for j, d in enumerate(cols):
        for i, k in enumerate(rows):
            if k in sig[d]: M[i, j] = sig[d][k][0]
    vmax = max(10.5, math.ceil(np.nanmax(np.abs(M))))
    cmap = LinearSegmentedColormap.from_list("gw", ["#4E9A51", "#FFFFFF", "#C0504D"]); cmap.set_bad("#EDEDED")
    axc = fig.add_axes([0.30, 0.215, 0.60, 0.46])
    im = axc.pcolormesh(np.ma.masked_invalid(M), cmap=cmap, vmin=-vmax, vmax=vmax, edgecolors="white", linewidth=1.5)
    axc.set_xlim(0, len(cols)); axc.set_ylim(len(rows), 0)
    axc.set_xticks(np.arange(len(cols))+.5); axc.set_xticklabels([lab(d) for d in cols], rotation=90, fontsize=13)
    axc.set_yticks(np.arange(len(rows))+.5)
    axc.set_yticklabels(["%s  (%s)" % (genes.get(k, "NA"), "novel" if cb.get(k, "novel") == "novel" else cb[k]) for k in rows], fontsize=13)
    axc.tick_params(length=0); axc.spines[:].set_visible(False)
    fig.text(0.115, 0.70, "(c)", fontsize=17, fontweight="bold")
    cax = fig.add_axes([0.915, 0.30, 0.012, 0.28]); cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("log$_2$FC (tumour vs. normal)", fontsize=13, rotation=90); cbar.ax.tick_params(labelsize=12); cbar.outline.set_linewidth(1.2)
    save(fig, "fig9")
    return len(rows), vmax, ge

fig7(); fig8(); r = fig9(); print("fig9c rows", r[0], "vmax", r[1], "tiers", r[2])
