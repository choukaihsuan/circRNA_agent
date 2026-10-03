"""Reproducible versions of thesis Figures 7, 8, 9(a,b,c) plus cross-dataset UpSet, for the 15-dataset main set
(figs/main_15ds) and the 19-dataset appendix set (figs/appendix_19ds).
Run: python3 make_figures.py   (needs matplotlib, numpy; inputs in ./data)
Layout/definitions follow the v19 captions: significant = pvalue<0.05 & |log2FC|>1 (edgeR_ciriquant).
"""
import os, collections, itertools, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
from matplotlib.colors import LinearSegmentedColormap
from thesis_common import *
summ, sig, tested, genes = load()
C_UP, C_DN, C_T1, C_T2 = "#d62728", "#2ca02c", "#1f77b4", "#ff7f0e"
def nsig(d): return len(sig[d])
def label(d): return "%s (%s)" % (d, META[d][1])
def tcount(d):
    c = collections.Counter(t for _, t in sig[d].values()); return c["Type_I"], c["Type_II"]

def fig7(ds, out, tag):
    order = sorted(ds, key=lambda d: -nsig(d)); n = len(order)
    fig = plt.figure(figsize=(11, 0.42*n+1.6)); ax = fig.add_axes([0.28, 0.12, 0.40, 0.80]); axp = fig.add_axes([0.70, 0.12, 0.28, 0.80])
    for i, d in enumerate(order):
        a, b = tcount(d); t = a+b; y = n-1-i
        ax.barh(y, 100*a/t, color=C_T1, edgecolor="white"); ax.barh(y, 100*b/t, left=100*a/t, color=C_T2, edgecolor="white")
        ax.text(101, y, "n=%d" % t, va="center", fontsize=7)
        r = 0.45*math.sqrt(t/max(nsig(x) for x in order))  # radius ∝ sqrt(count)
        axp.pie([a, b], colors=[C_T1, C_T2], radius=max(r, 0.06), center=(0.5, y), frame=True, wedgeprops=dict(edgecolor="white", linewidth=.3))
    ax.set_yticks(range(n)); ax.set_yticklabels([label(d) for d in order[::-1]], fontsize=8); ax.set_xlim(0, 112); ax.set_xlabel("% of significant DECs")
    fig.legend(handles=[plt.Rectangle((0,0),1,1,color=C_T1), plt.Rectangle((0,0),1,1,color=C_T2)], labels=["Type I", "Type II"], loc="lower center", ncol=2, fontsize=8, frameon=False)
    axp.set_xlim(0, 1); axp.set_ylim(-0.7, n-0.3); axp.axis("off"); axp.set_aspect("equal")
    fig.suptitle("Figure 7 – Type I/II distribution, %s" % tag, fontsize=10)
    fig.savefig(out+".png", dpi=200); fig.savefig(out+".pdf"); plt.close(fig)

def fig8(ds, out, tag):
    order = sorted(ds, key=lambda d: -nsig(d)); n = len(order)
    fig, ax = plt.subplots(figsize=(8, 0.42*n+1.4)); h = 0.38
    for i, d in enumerate(order):
        u = sum(1 for l, _ in sig[d].values() if l > 0); dn = nsig(d)-u; y = n-1-i
        if u: ax.barh(y+h/2, u, h, color=C_UP); ax.text(u*1.08, y+h/2, str(u), va="center", fontsize=7)
        if dn: ax.barh(y-h/2, dn, h, color=C_DN); ax.text(dn*1.08, y-h/2, str(dn), va="center", fontsize=7)
    ax.set_xscale("log"); ax.set_xlim(0.8, max(nsig(d) for d in order)*3)
    ax.set_yticks(range(n)); ax.set_yticklabels([label(d) for d in order[::-1]], fontsize=8)
    ax.set_xlabel("significant DECs (log10 scale; bars omitted when 0)"); ax.legend(handles=[plt.Rectangle((0,0),1,1,color=C_UP), plt.Rectangle((0,0),1,1,color=C_DN)], labels=["up in tumor", "down in tumor"], fontsize=8, loc="lower right")
    ax.set_title("Figure 8 – direction of dysregulation, %s" % tag, fontsize=10)
    fig.tight_layout(); fig.savefig(out+".png", dpi=200); fig.savefig(out+".pdf"); plt.close(fig)

def fig9(ds, out, tag):
    cnt = tiers(sig, ds); mx = max(cnt.values())
    dist = [sum(1 for v in cnt.values() if v == k) for k in range(1, mx+1)]
    fig, ax = plt.subplots(figsize=(4.5, 3.4)); ax.bar(range(1, mx+1), dist, color="#4c72b0"); ax.set_yscale("log")
    for k, v in enumerate(dist, 1): ax.text(k, v*1.1, str(v), ha="center", fontsize=8)
    ax.set_xlabel("number of datasets in which circRNA is significant"); ax.set_ylabel("circRNAs (log10)"); ax.set_title("Figure 9a – recurrence distribution, %s" % tag, fontsize=9)
    fig.tight_layout(); fig.savefig(out+"_9a.png", dpi=200); fig.savefig(out+"_9a.pdf"); plt.close(fig)
    # 9b summary table
    fig, ax = plt.subplots(figsize=(4.8, 2.2)); ax.axis("off")
    rows = [[">= %d datasets" % k, sum(1 for v in cnt.values() if v >= k)] for k in range(1, mx+1)]
    t = ax.table(cellText=rows, colLabels=["recurrence", "circRNAs"], loc="center"); t.scale(1, 1.3); ax.set_title("Figure 9b – recurrence summary, %s" % tag, fontsize=9)
    fig.savefig(out+"_9b.png", dpi=200, bbox_inches="tight"); fig.savefig(out+"_9b.pdf", bbox_inches="tight"); plt.close(fig)
    # 9c heatmap: circRNAs sig in >=4 datasets
    rowsid = sorted([k for k, v in cnt.items() if v >= 4], key=lambda k: (-cnt[k], k))
    cols = sorted(ds, key=lambda d: (META[d][0], META[d][1], d))
    M = np.full((len(rowsid), len(cols)), np.nan)
    for j, d in enumerate(cols):
        for i, k in enumerate(rowsid):
            if k in sig[d]: M[i, j] = sig[d][k][0]
    cmap = LinearSegmentedColormap.from_list("gwr", [C_DN, "white", C_UP]); cmap.set_bad("#cccccc")
    fig, ax = plt.subplots(figsize=(0.55*len(cols)+3.2, 0.3*len(rowsid)+2.6))
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-8, vmax=8, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(["%s · %s\n%s" % (META[d][0], META[d][1], d) for d in cols], rotation=90, fontsize=6.5)
    ax.set_yticks(range(len(rowsid))); ax.set_yticklabels(["%s  %s (%d)" % (genes.get(k, "NA"), k, cnt[k]) for k in rowsid], fontsize=6.5)
    fig.colorbar(im, ax=ax, label="log2FC (red up, green down; gray = not significant)", shrink=.6)
    ax.set_title("Figure 9c – circRNAs significant in ≥4 datasets (%d), %s" % (len(rowsid), tag), fontsize=9)
    fig.tight_layout(); fig.savefig(out+"_9c.png", dpi=200); fig.savefig(out+"_9c.pdf"); plt.close(fig)
    return rowsid

def upset(ds, out, tag, top=30):
    members = collections.defaultdict(set)
    for d in ds:
        for k in sig[d]: members[k].add(d)
    pat = collections.Counter(frozenset(v) for v in members.values() if len(v) >= 2)
    items = pat.most_common(top); order = sorted(ds, key=lambda d: -nsig(d)); n = len(order)
    fig = plt.figure(figsize=(0.32*top+4, 0.28*n+4))
    axb = fig.add_axes([0.30, 0.52, 0.68, 0.42]); axm = fig.add_axes([0.30, 0.06, 0.68, 0.44], sharex=axb); axs = fig.add_axes([0.03, 0.06, 0.12, 0.44])
    for x, (s, c) in enumerate(items):
        axb.bar(x, c, color="#1f77b4"); axb.text(x, c, str(c), ha="center", va="bottom", fontsize=6)
        ys = [n-1-order.index(d) for d in s]; axm.plot([x]*len(ys), ys, color="k", lw=1)
        for yy in range(n): axm.scatter(x, yy, s=14, color="k" if (n-1-yy) in [order.index(d) for d in s] else "#dddddd", zorder=3)
    axb.set_ylabel("circRNAs in exact intersection"); axb.set_xticks([]); axm.set_yticks(range(n)); axm.set_yticklabels([d for d in order[::-1]], fontsize=6.5); axm.set_xticks([])
    for sp in ("top", "right", "bottom"): axm.spines[sp].set_visible(False)
    axs.barh(range(n), [nsig(d) for d in order[::-1]], color="#888888"); axs.invert_xaxis(); axs.set_yticks([]); axs.set_xlabel("sig DECs per dataset", fontsize=7); axs.tick_params(labelsize=6)
    axb.set_title("UpSet – exact dataset intersections of circRNAs shared by ≥2 datasets (top %d), %s" % (top, tag), fontsize=9)
    fig.savefig(out+".png", dpi=200); fig.savefig(out+".pdf"); plt.close(fig)

for tag, ds, sub in (("15 dual-tool datasets (main)", DS15, "main_15ds"), ("all 19 datasets (appendix)", DS19, "appendix_19ds")):
    d = os.path.join(HERE, "figs", sub); os.makedirs(d, exist_ok=True)
    fig7(ds, os.path.join(d, "fig7_type_distribution"), tag); fig8(ds, os.path.join(d, "fig8_direction"), tag)
    rows = fig9(ds, os.path.join(d, "fig9"), tag); upset(ds, os.path.join(d, "upset"), tag)
    print(sub, "fig9c rows:", len(rows))
