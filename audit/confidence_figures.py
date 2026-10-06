"""Item 16 figures: Figure 10 (two panels) and the per-dataset report panels, from consensus_filter --dump-dist.

Reads existing dump TSVs only (<dump_root>/<GSE>/<dump_dirname>/<srr>.tsv); never runs consensus_filter and never
touches the pipeline. The same script and parameters are used for the draft and the final run; every run writes
<out_dir>/confidence_figures_params.json recording parameters, per-dataset counts and a hash of the inputs, so the
two runs can be compared.

A "locus" below is one row of a dump = one consensus circRNA in one sample (the same circRNA in two samples counts
twice). This is the unit behind the 628,793 figure.

usage: confidence_figures.py <dump_root> <out_dir> [--dump-dirname dump] [--seed 20261007] [--max-hexbin-points 0]
  --max-hexbin-points 0 (default) = no sampling, all loci are drawn.
Stops (exit 2) if any of the 15 datasets has a missing, empty or incomplete dump; nothing is substituted or dropped.
"""
import argparse
import base64
import glob
import hashlib
import io
import json
import os
import sys

import numpy as np
import pandas as pd

# expected number of per-sample dumps per dataset (analysed samples; GSE221107 = the 8 retained samples)
EXPECTED_SAMPLES = {
    "GSE108735": 14, "GSE113230": 6, "GSE121842": 6, "GSE130078": 12, "GSE133998": 12, "GSE136569": 10,
    "GSE143797": 8, "GSE148036": 10, "GSE171011": 8, "GSE221107": 8, "GSE248612": 12, "GSE77509": 12,
    "GSE97239": 6, "PRJNA553289": 12, "SRP156355": 12,
}
COLS = ["chrom", "start", "end", "n_supporting", "bsj_ciriquant", "dist_ciriquant", "bsj_dcc", "dist_dcc",
        "confidence_score"]
SLOP = 10  # consensus slop used for the run; agreement factor = 1 - dist/SLOP
AX_LO = 0.998  # panel (a) x-axis is truncated to [AX_LO, 1.0]; stated in the figure and caption
HIST_WIDTH = 0.25  # fixed confidence_score bin width, edges anchored at 0 (same bins for every dataset)


def load(root, dirname):
    problems, frames, digest = [], [], hashlib.sha256()
    for g, n_exp in EXPECTED_SAMPLES.items():
        files = sorted(f for f in glob.glob(os.path.join(root, g, dirname, "*.tsv")) if not f.endswith("_s.tsv"))
        if len(files) != n_exp:
            problems.append(f"{g}: {len(files)} dump files, expected {n_exp}")
        for f in files:
            try:
                d = pd.read_csv(f, sep="\t", na_values="NA")
            except Exception as e:  # unreadable / truncated
                problems.append(f"{g}/{os.path.basename(f)}: unreadable ({e})")
                continue
            if not set(COLS) <= set(d.columns):
                problems.append(f"{g}/{os.path.basename(f)}: missing columns {sorted(set(COLS) - set(d.columns))}")
                continue
            if d.confidence_score.notna().sum() == 0:
                problems.append(f"{g}/{os.path.basename(f)}: no consensus loci (empty or incomplete)")
                continue
            digest.update(f"{g}/{os.path.basename(f)}:{os.path.getsize(f)}\n".encode())
            d["dataset"] = g
            d["sample"] = os.path.splitext(os.path.basename(f))[0]
            frames.append(d)
    if problems:
        print("STOP: dump set is incomplete; nothing was plotted.", file=sys.stderr)
        for p in problems:
            print("  - " + p, file=sys.stderr)
        sys.exit(2)
    df = pd.concat(frames, ignore_index=True)
    df = df[df.confidence_score.notna()].copy()
    df["has_ciri"] = df.dist_ciriquant.notna()
    df["has_dcc"] = df.dist_dcc.notna()
    lc = np.log2(df.bsj_ciriquant.astype(float) + 1)
    ld = np.log2(df.bsj_dcc.astype(float) + 1)
    df["log_mean"] = np.where(df.has_ciri & df.has_dcc, (lc + ld) / 2, np.where(df.has_ciri, lc, ld))
    return df, digest.hexdigest()


def bucket_counts(d):
    d = d.dropna()
    return {"zero": int((d == 0).sum()), "one_three": int(((d >= 1) & (d <= 3)).sum()),
            "four_ten": int(((d >= 4) & (d <= 10)).sum()), "n": int(len(d))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump_root")
    ap.add_argument("out_dir")
    ap.add_argument("--dump-dirname", default="dump")
    ap.add_argument("--seed", type=int, default=20261007)
    ap.add_argument("--max-hexbin-points", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from scipy.stats import spearmanr

    df, in_hash = load(a.dump_root, a.dump_dirname)
    rng = np.random.default_rng(a.seed)

    # ── per-dataset numbers ──
    per = {}
    for g, sub in df.groupby("dataset"):
        b = bucket_counts(sub.loc[sub.has_dcc, "dist_dcc"])
        per[g] = {"loci": int(len(sub)), "dcc_supported": b["n"], "dist0": b["zero"],
                  "dist_gt0": b["n"] - b["zero"], "dist_1_3": b["one_three"], "dist_4_10": b["four_ten"],
                  "frac_dist0": b["zero"] / b["n"] if b["n"] else float("nan"),
                  "any_diff": int(((sub.dist_dcc > 0) | (sub.dist_ciriquant > 0)).sum())}
        per[g]["frac_identical"] = 1 - per[g]["any_diff"] / per[g]["loci"]
    order = sorted(per, key=lambda g: per[g]["loci"])  # ascending: smallest first, drawn bottom-up
    below = [g for g in order if per[g]["frac_identical"] < AX_LO]
    if below:
        print(f"STOP: dist_dcc=0 fraction < {AX_LO} in {below}; the truncated axis would hide them. "
              f"Widen AX_LO deliberately and change the caption.", file=sys.stderr)
        sys.exit(2)

    # ── Figure 10 ──
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(14, 6.4), gridspec_kw={"width_ratios": [1.05, 1]})
    ys = np.arange(len(order))
    xs = [per[g]["frac_identical"] for g in order]
    axa.plot(xs, ys, "o", color="#4C72B0", ms=7)
    for y, g in zip(ys, order):
        p = per[g]
        axa.annotate(f"{p['any_diff']:,} of {p['loci']:,}", (p["frac_identical"], y),
                     xytext=(8, 0), textcoords="offset points", va="center", fontsize=8)
    axa.set_yticks(ys)
    axa.set_yticklabels(order, fontsize=9)
    axa.set_xlim(AX_LO, 1.0)
    axa.set_xlabel(f"Fraction of loci with identical coordinates from both callers\n[x-axis truncated to {AX_LO:.3f}\u20131.000]")
    axa.set_title("(a) Coordinate agreement per dataset (sorted by locus count)", fontsize=10)
    axa.grid(axis="x", color="#e5e7eb", lw=0.6)
    axa.set_axisbelow(True)
    axa.text(0.0, -0.19, "label = loci with different coordinates (n of all loci in the dataset)", transform=axa.transAxes,
             fontsize=8, color="#555")

    pts = df
    sampled_note = "all loci drawn, no sampling"
    if a.max_hexbin_points and len(df) > a.max_hexbin_points:
        pts = df.iloc[np.sort(rng.choice(len(df), a.max_hexbin_points, replace=False))]
        sampled_note = f"random sample {len(pts):,} of {len(df):,} (seed {a.seed})"
    hb = axb.hexbin(pts.log_mean, pts.confidence_score, gridsize=70, cmap="viridis", mincnt=1,
                    norm=LogNorm())
    lim = [0, max(pts.log_mean.max(), pts.confidence_score.max()) * 1.02]
    axb.plot(lim, lim, "--", color="#d62728", lw=1.2, label="y = x")
    rho_all = spearmanr(df.confidence_score, df.log_mean).correlation
    two = df[df.has_ciri & df.has_dcc]
    rho_two = spearmanr(two.confidence_score, two.log_mean).correlation
    axb.text(0.04, 0.95, f"\u03c1 = {rho_two:.4f}\nn = {len(df):,} consensus loci",
             transform=axb.transAxes, va="top", fontsize=9, bbox=dict(fc="white", ec="#bbb", alpha=0.9))
    axb.set_xlabel("mean log2(BSJ + 1) over supporting tools")
    axb.set_ylabel("confidence_score")
    axb.set_title("(b) confidence_score vs mean log2(BSJ + 1)", fontsize=10)
    axb.legend(loc="lower right", fontsize=8)
    fig.colorbar(hb, ax=axb, label="loci per hexagon (log scale)")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    f10 = os.path.join(a.out_dir, "figure10.png")
    fig.savefig(f10, dpi=200)
    fig.savefig(os.path.join(a.out_dir, "figure10.pdf"))
    plt.close(fig)

    # ── per-dataset report panels ──
    rep_dir = os.path.join(a.out_dir, "report_panels")
    os.makedirs(rep_dir, exist_ok=True)
    for g, sub in df.groupby("dataset"):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.2), gridspec_kw={"width_ratios": [1.6, 1]})
        a1.hist(sub.confidence_score, bins=np.arange(np.floor(sub.confidence_score.min() / HIST_WIDTH) * HIST_WIDTH, np.ceil(sub.confidence_score.max() / HIST_WIDTH) * HIST_WIDTH + HIST_WIDTH, HIST_WIDTH), color="#4C72B0")
        a1.set_xlabel("Confidence score")
        a1.set_ylabel("Loci")
        a1.set_title("Confidence score of each circRNA call", fontsize=10)
        p = per[g]
        vals = [p["dist0"], p["dist_1_3"], p["dist_4_10"]]
        bars = a2.bar(["0", "1–3", "4–10"], vals, color=["#4C72B0", "#DD8452", "#C44E52"])
        for bar, v in zip(bars, vals):
            a2.annotate(f"{v:,}\n({v / p['dcc_supported']:.2%})", (bar.get_x() + bar.get_width() / 2, v),
                        ha="center", va="bottom", fontsize=8)
        a2.set_ylim(0, p["dcc_supported"] * 1.28)
        a2.set_xlabel("Coordinate difference between callers (bp)")
        a2.set_ylabel("Loci")
        a2.set_title("Do the callers agree?", fontsize=10)
        fig.suptitle(f"{g}: {p['loci']:,} loci", fontsize=10)
        fig.tight_layout()
        png = os.path.join(rep_dir, f"{g}_confidence.png")
        fig.savefig(png, dpi=130)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=130)
        plt.close(fig)
        b64 = base64.b64encode(buf.getvalue()).decode()
        note_en = ("Left: distribution of confidence_score. Right: coordinate distance between the DCC call and the "
                   "representative coordinate. The score is shown together with the distance because almost all loci "
                   "have distance 0: the agreement term is near-constant, so the score varies with the BSJ read "
                   "term only. The 0 bar means the DCC coordinate equals the representative coordinate; it also "
                   "includes the rare loci where the CIRIquant coordinate differs instead.")
        note_zh = ("左：confidence_score 分布。右：DCC 座標與代表座標的距離。兩者並列，是因為幾乎所有 locus 的距離為 0，"
                   "一致性項近乎常數，分數的變化實際上只來自 BSJ reads 項。距離 0 的柱表示 DCC 座標與代表座標相同；"
                   "也包含極少數改為 CIRIquant 座標不同的 locus（見論文資料集層級說明）。")
        html = (f'<div class="conf-structure"><h3 data-en="Consensus confidence score and coordinate agreement">'
                f'共識信心分數與座標一致性</h3><img alt="{g} confidence score" style="max-width:100%" '
                f'src="data:image/png;base64,{b64}"><p data-en="{note_en}">{note_zh}</p></div>\n')
        open(os.path.join(rep_dir, f"{g}_confidence.html"), "w", encoding="utf-8").write(html)

    n_gt0 = sum(p["dist_gt0"] for p in per.values())
    n_ciri_only = int(((df.dist_dcc == 0) & (df.dist_ciriquant > 0)).sum())  # DCC coordinate was the representative
    n_any = int(((df.dist_dcc > 0) | (df.dist_ciriquant > 0)).sum())
    assert len(two) == len(df), "caption assumes every locus is supported by both callers"
    lo, hi = min(order, key=lambda g: per[g]["frac_identical"]), max(order, key=lambda g: per[g]["frac_identical"])
    cap = (f"Figure 10. (a) Per dataset, fraction of consensus loci at which both callers report identical coordinates: {per[lo]['frac_identical']:.4f} in {lo} to {per[hi]['frac_identical']:.4f} in {hi}. "
           f"The x-axis is truncated to {AX_LO:.3f}\u20131.000 (full range 0\u20131). Labels give the number of loci with different coordinates "
           f"out of all loci in that dataset. Datasets are sorted by locus count. A locus is one consensus circRNA in one sample. "
           f"(b) confidence_score against the mean log2(BSJ+1) over the two callers for all {len(df):,} loci (hexbin, log colour "
           f"scale; {sampled_note}); dashed line y = x; Spearman rho = {rho_two:.4f}. Every locus is supported by both callers. "
           f"The small cluster below the identity line comprises the {n_any:,} loci ({n_any / len(df):.2%}) at which the two "
           f"callers reported different coordinates ({n_gt0:,} with dist_dcc > 0, plus {n_ciri_only:,} with dist_ciriquant > 0 "
           f"where the DCC coordinate was the representative); these account for all of the visible departure from y = x.\n")
    open(os.path.join(a.out_dir, "figure10_caption.md"), "w", encoding="utf-8").write(cap)

    meta = {"dump_root": os.path.abspath(a.dump_root), "dump_dirname": a.dump_dirname, "seed": a.seed,
            "max_hexbin_points": a.max_hexbin_points, "hexbin_sampling": sampled_note, "hexbin_gridsize": 70,
            "slop": SLOP, "axis_a_xlim": [AX_LO, 1.0], "total_loci": int(len(df)), "loci_any_dist_gt0": n_any, "loci_dist_dcc_gt0": n_gt0, "loci_dist_ciri_gt0_only": n_ciri_only,
            "rho_all": float(rho_all), "rho_two_tool": float(rho_two), "input_hash_sha256": in_hash,
            "per_dataset": per, "dataset_order_bottom_to_top": order}
    json.dump(meta, open(os.path.join(a.out_dir, "confidence_figures_params.json"), "w"), indent=1)
    print(f"loci {len(df):,}  rho_two {rho_two:.4f}  input_hash {in_hash[:12]}  {sampled_note}")
    for g in reversed(order):
        p = per[g]
        print(f"{g:12s} loci {p['loci']:>8,}  dcc {p['dcc_supported']:>8,}  dist_dcc>0 {p['dist_gt0']:>4,}  any {p['any_diff']:>4,}  "
              f"ident {p['frac_identical']:.5f}")


if __name__ == "__main__":
    main()
