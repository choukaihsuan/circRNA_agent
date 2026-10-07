"""Figure 12: composite biomarker score distribution (GSE133998, n=79 significant circRNAs).
Reads fig12_gse133998_components.tsv (from fig12_components.py) and fig12_medians.json.
usage: make_fig12.py <HOME> <out_dir>
"""
import json
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

HOME, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(os.path.join(HOME, "item16", "fig12_gse133998_components.tsv"), sep="\t")
medians = json.load(open(os.path.join(HOME, "item16", "fig12_medians.json")))

COMP = ["sig", "fc", "conf", "known", "mirna", "rbp"]
LETTER = {"sig": "S", "fc": "F", "conf": "C", "known": "K", "mirna": "M", "rbp": "R"}

fig = plt.figure(figsize=(5.5, 5.5), constrained_layout=True)
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.08, 0.92])

# (a) composite score histogram
axa = fig.add_subplot(gs[0])
axa.hist(df["biomarker_score"], bins=20, color="#4C72B0", edgecolor="white", lw=0.4)
axa.set_xlabel("Composite biomarker_score (equal-weight default)")
axa.set_ylabel("circRNAs")
axa.set_title("(a)", loc="left", fontweight="bold")
axa.text(0.98, 0.95, f"n = {len(df):,}", transform=axa.transAxes, ha="right", va="top",
         bbox=dict(fc="white", ec="#bbb", alpha=0.9))

# label row for (b), separate from the subplot titles below it so nothing is clipped or collides
axlbl = fig.add_subplot(gs[1])
axlbl.axis("off")
axlbl.text(0.0, 0.5, "(b)", ha="left", va="center", fontweight="bold", transform=axlbl.transAxes)

# (b) six components, 1x6 small multiples; K drawn as a 0/1 bar chart, not a histogram
gs_b = gs[2].subgridspec(1, 6, wspace=0.15)
for i, c in enumerate(COMP):
    ax = fig.add_subplot(gs_b[i])
    vals = df[c].dropna()
    if c == "known":
        counts = [int((vals == 0).sum()), int((vals == 1).sum())]
        ax.bar(["0", "1"], counts, color="#C44E52", width=0.6)
        ax.set_ylim(0, max(counts) * 1.15)
        ax.set_title("K (binary)", fontsize=8)
        for x, v in enumerate(counts):
            ax.annotate(str(v), (x, v), ha="center", va="bottom", fontsize=7)
    else:
        ax.hist(vals, bins=np.linspace(0, 1, 11), color="#4C72B0", orientation="vertical")
        ax.set_xlim(0, 1)
        ax.set_title(LETTER[c], fontsize=8)
    if i == 0:
        ax.set_ylabel("circRNAs", fontsize=7)
    ax.tick_params(labelsize=6)
    if c != "known":
        ax.set_xticks([0, 1])

fig.savefig(os.path.join(OUT, "figure12.png"), dpi=600)
fig.savefig(os.path.join(OUT, "figure12.pdf"), dpi=600)
plt.close(fig)

med_tbl = []
for g, m in medians.items():
    med_tbl.append((g, m["n_sig"], m["median"], m["top"]))
top_full = [LETTER[m["top"]] for m in medians.values()]
from collections import Counter
tally_full = Counter(top_full)
top_exclK = []
for g, m in medians.items():
    med = m["median"]
    c5 = {k: v for k, v in med.items() if k != "known"}
    top_exclK.append(LETTER[max(c5, key=c5.get)])
tally_exclK = Counter(top_exclK)
min_tally = Counter()
for g, m in medians.items():
    med = m["median"]
    c5 = {k: v for k, v in med.items() if k != "known"}
    min_tally[LETTER[min(c5, key=c5.get)]] += 1

cap = f"""Figure 12. Composite biomarker_score and its six components (GSE133998, n = {len(df)} significant circRNAs, the same dataset as Figure 11 / Table 5 / §3.2).

(a) Distribution of the composite score (equal-weight default, §2.8).

(b) The six components: S (sig_norm), F (fc_norm), C (conf_norm), M (mirna_norm), R (rbp_norm) are shown as histograms of their within-dataset min–max normalized value (0–1, 10 bins). K (known_bonus) is a raw, unnormalized 0/1 flag (in circBase or not) and is shown as a 2-bar count, not a histogram, to avoid implying it is continuous.

Scope: all values are normalized within GSE133998's own significant set (§2.8); they are not comparable across datasets and this figure is not a 15-dataset pooled histogram.

Ceiling effects: S (capped at −log10 p = 10) and F (capped at |log2FC| = 5) are clipped before normalization, so circRNAs at or above the cap pile up at 1.0; this spike is a design artifact of the cap, not a biological signal.

K is unnormalized and binary (0/1); its "distribution" is two bars, not a continuous shape.

Representativeness across the 15 analysed datasets (median of each normalized component in the dataset's own significant set): by raw median, K has the highest value in {tally_full.get('K',0)}/15 datasets (including GSE133998) -- but this is a structural property of a 0/1 variable with >50% of calls in circBase, not evidence that K is more discriminative than the continuous components. Excluding K, F has the highest median in {tally_exclK.get('F',0)}/15 datasets (GSE133998 included), R in {tally_exclK.get('R',0)}/15, C in {tally_exclK.get('C',0)}/15 -- GSE133998's own ranking matches the modal pattern but is not universal. M (mirna_norm) has the lowest median of the five continuous components in {min_tally.get('M',0)}/15 datasets, consistent with §4.2/§4.10's claim that the interaction-count components carry limited rank information for most datasets.

Equal weights are the default (§2.8); weights are configurable per component.
"""
open(os.path.join(OUT, "figure12_caption.md"), "w", encoding="utf-8").write(cap)
print("wrote figure12.png/pdf/caption.md to", OUT)
