"""Reproducible versions of thesis Figures 7, 8, 9(a,b,c) plus cross-dataset UpSet, for the 15-dataset main set
(figs/main_15ds) and the 19-dataset appendix set (figs/appendix_19ds). No in-figure titles (titles belong in captions).
Run: python3 make_figures.py   (needs matplotlib, numpy; inputs in ./data)
Significant = pvalue<0.05 & |log2FC|>1 (edgeR_ciriquant), as in thesis Sec 2.7.
"""
import os, collections, math
import matplotlib; matplotlib.use("Agg"); import matplotlib.ticker
import matplotlib.pyplot as plt, numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle, Patch
from thesis_common import *
plt.rcParams.update({"font.family": "DejaVu Serif", "pdf.fonttype": 42, "hatch.linewidth": 0.8})
summ, sig, tested, genes = load()
C_UP, C_DN, C_T1, C_T2 = "#d62728", "#2ca02c", "#1f77b4", "#ff7f0e"
SMALL_N = 20
def nsig(d): return len(sig[d])
def label(d): return "%s (%s)" % (d, META[d][1])
def tcount(d):
    c = collections.Counter(t for _, t in sig[d].values()); return c["Type_I"], c["Type_II"]
def save(fig, out, **kw):
    fig.savefig(out+".png", dpi=200, **kw); fig.savefig(out+".pdf", **kw); plt.close(fig)

def fig7(ds, out):
    order = sorted(ds, key=lambda d: -nsig(d)); n = len(order)
    fig, ax = plt.subplots(figsize=(8.2, 0.40*n+1.6)); ax.set_position([0.36, 0.14, 0.52, 0.82])
    for i, d in enumerate(order):
        a, b = tcount(d); t = a+b; y = n-1-i; small = t < SMALL_N; hp = "///" if small else None
        ax.barh(y, 100*a/t, color=C_T1, edgecolor="white", hatch=hp); ax.barh(y, 100*b/t, left=100*a/t, color=C_T2, edgecolor="white", hatch=hp)
        ax.text(102, y, "n=%d%s" % (t, " †" if small else ""), va="center", fontsize=8, fontweight="bold" if small else "normal")
    ax.set_yticks(range(n)); ax.set_yticklabels([label(d) for d in order[::-1]], fontsize=8.5); ax.set_xlim(0, 100); ax.set_ylim(-0.6, n-0.4)
    ax.set_xlabel("% of significant DECs"); ax.spines[["top", "right"]].set_visible(False)
    hs = [Patch(color=C_T1), Patch(color=C_T2), Patch(facecolor="#bbbbbb", edgecolor="white", hatch="///")]
    fig.legend(hs, ["Type I", "Type II", "† n < %d (percentages unstable)" % SMALL_N], loc="lower center", ncol=3, fontsize=8, frameon=False)
    save(fig, out)

def fig8(ds, out):
    def up(d): return sum(1 for l, _ in sig[d].values() if l > 0)
    order = sorted(ds, key=lambda d: (-up(d)/nsig(d), -nsig(d))); n = len(order)
    fig, ax = plt.subplots(figsize=(8.4, 0.40*n+1.6)); ax.set_position([0.36, 0.17, 0.46, 0.80])
    for i, d in enumerate(order):
        u = up(d); dn = nsig(d)-u; t = nsig(d); y = n-1-i; pu = 100*u/t; small_note = []
        ax.barh(y, pu, color=C_UP, edgecolor="white"); ax.barh(y, 100-pu, left=pu, color=C_DN, edgecolor="white")
        for val, x0, x1, col in ((u, 0, pu, C_UP), (dn, pu, 100, C_DN)):
            if val == 0:  # explicit zero marker at that end
                xm = 2.5 if col == C_UP else 97.5
                ax.text(xm, y, "0", ha="center", va="center", fontsize=8, color=col, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=col, lw=1))
            elif x1-x0 >= 7: ax.text((x0+x1)/2, y, str(val), ha="center", va="center", fontsize=7.5, color="white")
            else: small_note.append("%s%d" % ("↑" if col == C_UP else "↓", val))   # too narrow to label inside the bar
        ax.text(101.5, y, "n=%d%s%s" % (t, " †" if t < SMALL_N else "", " (%s)" % small_note[0] if small_note else ""), va="center", fontsize=8)
    ax.axvline(50, color="k", lw=0.6, ls=":")
    ax.set_yticks(range(n)); ax.set_yticklabels([label(d) for d in order[::-1]], fontsize=8.5); ax.set_xlim(0, 100); ax.set_ylim(-0.6, n-0.4)
    ax.set_xlabel("% of significant DECs"); ax.spines[["top", "right"]].set_visible(False)
    fig.legend([Patch(color=C_UP), Patch(color=C_DN)], ["up in tumor", "down in tumor"], loc="lower center", bbox_to_anchor=(0.5, 0.035), ncol=2, fontsize=8, frameon=False)
    fig.text(0.36, 0.008, "Sorted by share up-regulated; boxed 0 = no DECs in that direction; † n < %d" % SMALL_N, fontsize=7)
    save(fig, out)

def fig9(ds, out):
    cnt = tiers(sig, ds); mx = max(cnt.values())
    dist = [sum(1 for v in cnt.values() if v == k) for k in range(1, mx+1)]
    fig, ax = plt.subplots(figsize=(4.5, 3.4)); ax.bar(range(1, mx+1), dist, color="#4c72b0"); ax.set_yscale("log")
    for k, v in enumerate(dist, 1): ax.text(k, v*1.1, str(v), ha="center", fontsize=8)
    ax.set_xlabel("number of datasets in which circRNA is significant"); ax.set_ylabel("circRNAs (log10)"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); save(fig, out+"_9a")
    fig, ax = plt.subplots(figsize=(3.6, 0.35*mx+0.9)); ax.axis("off")
    rows = [[">= %d datasets" % k, sum(1 for v in cnt.values() if v >= k)] for k in range(1, mx+1)]
    t = ax.table(cellText=rows, colLabels=["recurrence", "circRNAs"], loc="center"); t.scale(1, 1.3)
    save(fig, out+"_9b", bbox_inches="tight")
    rowsid = sorted([k for k, v in cnt.items() if v >= 4], key=lambda k: (-cnt[k], k))
    cols = sorted(ds, key=lambda d: (META[d][0], META[d][1], d))
    M = np.full((len(rowsid), len(cols)), np.nan)
    for j, d in enumerate(cols):
        for i, k in enumerate(rowsid):
            if k in sig[d]: M[i, j] = sig[d][k][0]
    cmap = LinearSegmentedColormap.from_list("gwr", [C_DN, "white", C_UP])
    fig, ax = plt.subplots(figsize=(0.55*len(cols)+3.6, 0.3*len(rowsid)+2.6))
    for i in range(len(rowsid)):   # hatched dark-gray cells = not significant (distinct from pale mid-scale colors)
        for j in range(len(cols)):
            if np.isnan(M[i, j]):
                ax.add_patch(Rectangle((j-.5, i-.5), 1, 1, facecolor="#9e9e9e", edgecolor="#4d4d4d", hatch="////", lw=0))
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-8, vmax=8, aspect="auto")
    ax.set_xlim(-.5, len(cols)-.5); ax.set_ylim(len(rowsid)-.5, -.5)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(["%s · %s\n%s" % (META[d][0], META[d][1], d) for d in cols], rotation=90, fontsize=7)
    ax.set_yticks(range(len(rowsid))); ax.set_yticklabels(["%s  %s (%d)" % (genes.get(k, "NA"), k, cnt[k]) for k in rowsid], fontsize=7)
    cb = fig.colorbar(im, ax=ax, label="log2FC (red up, green down)", shrink=.6)
    fig.legend([Patch(facecolor="#9e9e9e", edgecolor="#4d4d4d", hatch="////")], ["not significant"], loc="lower right", fontsize=8, frameon=False)
    fig.tight_layout(); save(fig, out+"_9c")
    return rowsid

def upset(ds, out, top=30):
    members = collections.defaultdict(set)
    for d in ds:
        for k in sig[d]: members[k].add(d)
    pat = collections.Counter(frozenset(v) for v in members.values() if len(v) >= 2)
    items = sorted(pat.items(), key=lambda x: (-len(x[0]), -x[1], sorted(x[0])))[:top]   # by degree, then size
    order = sorted(ds, key=lambda d: -nsig(d)); n = len(order); pos = {d: n-1-i for i, d in enumerate(order)}
    fig = plt.figure(figsize=(0.32*top+4, 0.28*n+4.2))
    axb = fig.add_axes([0.30, 0.54, 0.68, 0.42]); axm = fig.add_axes([0.30, 0.10, 0.68, 0.42], sharex=axb); axs = fig.add_axes([0.03, 0.10, 0.12, 0.42])
    for x, (s, c) in enumerate(items):
        axb.bar(x, c, color="#1f77b4"); axb.text(x, c, str(c), ha="center", va="bottom", fontsize=6.5)
        ys = [pos[d] for d in s]; axm.plot([x]*2, [min(ys), max(ys)], color="k", lw=1)
        axm.scatter([x]*n, range(n), s=14, color=["k" if yy in ys else "#dddddd" for yy in range(n)], zorder=3)
    axb.set_ylabel("circRNAs in exact intersection"); axb.tick_params(labelbottom=False); axb.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True)); axb.spines[["top", "right"]].set_visible(False)
    axm.set_yticks(range(n)); axm.set_yticklabels(order[::-1], fontsize=7); axm.set_xticks(range(len(items))); axm.set_xticklabels([str(len(s)) for s, _ in items], fontsize=7)
    axm.set_xlabel("number of datasets in intersection (sorted by degree, then size)", fontsize=8)
    for sp in ("top", "right", "bottom"): axm.spines[sp].set_visible(False)
    axs.barh(range(n), [nsig(d) for d in order[::-1]], color="#888888"); axs.invert_xaxis(); axs.set_yticks([]); axs.set_xlabel("sig DECs per dataset", fontsize=7); axs.tick_params(labelsize=6); axs.spines[["top", "left"]].set_visible(False)
    save(fig, out)

for ds, sub in ((DS15, "main_15ds"), (DS19, "appendix_19ds")):
    d = os.path.join(HERE, "figs", sub); os.makedirs(d, exist_ok=True)
    fig7(ds, os.path.join(d, "fig7_type_distribution")); fig8(ds, os.path.join(d, "fig8_direction"))
    rows = fig9(ds, os.path.join(d, "fig9")); upset(ds, os.path.join(d, "upset"))
    print(sub, "fig9c rows:", len(rows))
