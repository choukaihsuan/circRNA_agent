"""Item 16A Figure 13 (standalone), stacked layout: (a) ROC curves on top, (b) AUROC/CI forest plot
on the bottom. Third layout pass -- same data, same AUROC/CI/DeLong numbers, same curves; only the
panel arrangement, panel letters, legend placement and font scale change.

Reads only ~/item16a/joined_scores.tsv (already-built join of SRR444655 --dump-dist output to the
live SRP011042 ground truth) and the DeLong AUROC/CI numbers already computed in R/pROC (hardcoded
below from ablation_auroc.R's stdout, verbatim -- no new AUROC/DeLong computation here).

Output: 3300 x 3800 px at 300 dpi -> figsize (11.0, 12.6667) in, exactly what the pixel/dpi target
forces (there is no freedom to choose a narrower native width while keeping both the pixel width
and the dpi fixed).

Font-size note (worth flagging explicitly): this figure prints at a FIXED 5.9-inch page width, so
the native-to-print scale factor is 5.9/11.0 = 0.5364 regardless of how the native HEIGHT is chosen
-- changing only the figsize height (as the two previous passes did, going from figsize height 4.6in
to 7in while keeping width at 11in and rcParams font.size=8) does not change that width-driven scale
factor at all, so it could not by itself have pushed effective print text past ~8*0.5364=4.3pt. To
actually hit ">=9pt at 5.9-inch print width" given native width is pinned at 11.0in, the rcParams
point size itself has to go up by about 11/5.9=1.864x: 9pt / 0.5364 = 16.78pt needed at minimum, so
this version uses a base rcParams font.size of 18pt (prints at 18*0.5364=9.65pt, a small margin
above the 9pt floor). Marker/line widths are left as in the previous pass (B6: styles unchanged).

usage: python3 make_figure13_stacked.py <joined_scores.tsv> <out_png>
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

BASE_PT = 18  # -> 18 * (5.9/11.0) = 9.65pt at 5.9-inch print width (>=9pt target, see module docstring)
matplotlib.rcParams.update({
    "font.size": BASE_PT, "axes.titlesize": BASE_PT, "axes.labelsize": BASE_PT,
    "xtick.labelsize": BASE_PT, "ytick.labelsize": BASE_PT, "legend.fontsize": BASE_PT,
})
import matplotlib.pyplot as plt

IN_TSV, OUT_PNG = sys.argv[1], sys.argv[2]
FIG_W_IN, FIG_H_IN = 3300 / 300, 3800 / 300  # (11.0, 12.6667) -- forced by the 3300x3800@300dpi spec
DPI = 300

# AUROC + 95% CI (DeLong, via pROC in the isolated `rauc` env), label version (i) ER>=1.5/<=0.5,
# n_labeled=621 (558 pos / 63 neg) -- from ablation_auroc.R stdout, reproduced verbatim here.
RESULTS_I = {
    "S1": {"label": "S1: confidence_score (full)",   "auc": 0.6285, "lo": 0.5519, "hi": 0.7052, "const": False},
    "S2": {"label": "S2: mean log2(BSJ+1)",           "auc": 0.6260, "lo": 0.5477, "hi": 0.7043, "const": False},
    "S3": {"label": "S3: mean(1 - dist/slop)",        "auc": 0.5044, "lo": 0.4884, "hi": 0.5205, "const": False},
    "S4": {"label": "S4: n_supporting (constant)",    "auc": 0.5000, "lo": None,   "hi": None,   "const": True},
}
DELONG_S1_S2_P = 0.551992
DELONG_S1_S3_P = 0.00240244

STYLE = {
    "S1": dict(color="#4C72B0", ls="-",  marker="o", lw=1.6, ms=4),
    "S2": dict(color="#DD8452", ls="--", marker="s", lw=1.6, ms=4),
    "S3": dict(color="#C44E52", ls="-.", marker="^", lw=1.6, ms=4),
    "S4": dict(color="#888888", ls=":",  marker="D", lw=1.2, ms=4),
}


def roc_points(score: np.ndarray, label: np.ndarray):
    """Empirical ROC curve (TPR vs FPR), higher score = more likely positive (direction '<' in pROC
    terms: controls have lower values). Threshold-sweep over every observed score value, ties handled
    by evaluating thresholds at each unique score (step function, consistent with a rank-based AUC)."""
    order = np.argsort(-score, kind="mergesort")
    s, y = score[order], label[order]
    P, N = int((y == 1).sum()), int((y == 0).sum())
    tp = np.cumsum(y == 1)
    fp = np.cumsum(y == 0)
    uniq_end = np.r_[np.where(np.diff(s) != 0)[0], len(s) - 1]
    tpr = np.r_[0.0, tp[uniq_end] / P]
    fpr = np.r_[0.0, fp[uniq_end] / N]
    return fpr, tpr


def main():
    d = pd.read_csv(IN_TSV, sep="\t")
    sub = d[d["label_i"] != -1].copy()
    y = sub["label_i"].values.astype(int)
    assert len(sub) == 621 and int((y == 1).sum()) == 558 and int((y == 0).sum()) == 63, \
        f"unexpected n for label (i): n={len(sub)} pos={int((y==1).sum())} neg={int((y==0).sum())}"

    fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN), constrained_layout=True)
    # top panel taller / closer to square: at 18pt base font a 5-line legend needs real room, and a
    # squat wide ROC panel leaves it nowhere to go without covering the curves (first attempt at
    # 1.05:0.95 did exactly that -- legend swallowed the whole lower-left of the curve area).
    gs = fig.add_gridspec(2, 1, height_ratios=[1.35, 0.65])
    axa = fig.add_subplot(gs[0])  # (a) ROC curves
    axb = fig.add_subplot(gs[1])  # (b) forest plot

    # ── (a) ROC curves ──
    # Legend text shortened (AUROC values dropped from these labels): at 18pt base font, the full
    # "Sn: <formula> (AUROC=X.XXX)" strings from the previous (side-by-side) layout made the legend
    # box wider/taller than the plot area itself. The AUROC values are not lost -- panel (b) is the
    # authoritative display of every AUROC + its CI; the legend here only needs to say which curve is
    # which.
    axa.plot([0, 1], [0, 1], "--", color="#999999", lw=1.0, label="chance")
    SHORT_LABEL = {"S1": "S1: confidence_score (full)", "S2": "S2: mean log2(BSJ+1)",
                   "S3": "S3: mean(1 − dist/slop)"}
    for nm, col in [("S1", "S1_confidence_score"), ("S2", "S2_mean_log2bsj"), ("S3", "S3_mean_agreement")]:
        fpr, tpr = roc_points(sub[col].values.astype(float), y)
        st = STYLE[nm]
        axa.plot(fpr, tpr, color=st["color"], ls=st["ls"], lw=st["lw"], label=SHORT_LABEL[nm])
    # S4: constant predictor -> ROC curve degenerates to the diagonal; drawn as a distinct dotted
    # overlay directly on the chance line (unchanged from previous passes)
    st4 = STYLE["S4"]
    axa.plot([0, 1], [0, 1], ls=st4["ls"], color=st4["color"], lw=st4["lw"],
              label="S4: n_supporting (constant)")
    axa.set_xlim(-0.02, 1.02)
    axa.set_ylim(-0.02, 1.02)
    axa.set_xlabel("False positive rate")
    axa.set_ylabel("True positive rate")
    axa.set_title("(a)", loc="left", fontweight="bold")
    axa.grid(color="#e5e7eb", lw=0.6)
    axa.set_axisbelow(True)
    axa.text(0.03, 0.97, "label (i): ER≥1.5/≤0.5\nn=621 (558 pos / 63 neg)\nSRR444655 (SRP011042)",
              transform=axa.transAxes, va="top", ha="left",
              bbox=dict(fc="white", ec="#bbb", alpha=0.9))
    # B5: legend back inside the ROC panel, lower-right corner (stacked layout gives it room again)
    axa.legend(loc="lower right", framealpha=0.9)

    # ── (b) AUROC + 95% CI forest plot ──
    order_names = ["S1", "S2", "S3", "S4"]
    ys = np.arange(len(order_names))[::-1]  # S1 at top
    axb.axvline(0.5, ls="--", color="#999999", lw=1.0, zorder=0)
    for y_pos, nm in zip(ys, order_names):
        r = RESULTS_I[nm]
        st = STYLE[nm]
        if r["const"]:
            axb.plot([r["auc"]], [y_pos], marker=st["marker"], mfc="white", mec=st["color"],
                      mew=1.6, ms=9, ls="none", zorder=3)
        else:
            axb.plot([r["lo"], r["hi"]], [y_pos, y_pos], color=st["color"], lw=2.0, zorder=2)
            axb.plot([r["auc"]], [y_pos], marker=st["marker"], mfc=st["color"], mec=st["color"],
                      ms=8, ls="none", zorder=3)
    axb.set_yticks(ys)
    axb.set_yticklabels([RESULTS_I[nm]["label"] for nm in order_names])
    axb.set_xlim(0.40, 0.90)
    axb.set_xlabel("AUROC (DeLong 95% CI where defined)")
    axb.set_title("(b)", loc="left", fontweight="bold")
    axb.grid(axis="x", color="#e5e7eb", lw=0.6)
    axb.set_axisbelow(True)

    # DeLong p-value annotation between the S1 and S2 rows (position/content unchanged)
    y_s1, y_s2 = ys[0], ys[1]
    y_mid = (y_s1 + y_s2) / 2
    axb.annotate("", xy=(0.80, y_s1), xytext=(0.80, y_s2),
                 arrowprops=dict(arrowstyle="-", color="#555555", lw=0.8), annotation_clip=False)
    axb.text(0.815, y_mid, f"DeLong\np={DELONG_S1_S2_P:.3f}", va="center", ha="left", color="#333333")

    fig.savefig(OUT_PNG, dpi=DPI)
    plt.close(fig)

    from PIL import Image
    im = Image.open(OUT_PNG)
    print(f"saved {OUT_PNG}  size={im.size}  dpi={im.info.get('dpi')}  "
          f"print_width_in=5.9  print_height_in={5.9 * im.size[1] / im.size[0]:.3f}")


if __name__ == "__main__":
    main()
