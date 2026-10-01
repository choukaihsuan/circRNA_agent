"""
run_permutation.py -- Set permutation test for cross-dataset circRNA recurrence.

Question: under the null hypothesis of random, independent per-dataset
significance calls, how many circRNAs would we expect to be "significant"
in >=2 / >=3 / >=4 / >=5 datasets just by chance? Is the observed
530 / 118 / 31 / 7 (out of 19 tumor-vs-normal datasets, GSE323364 cell-line
excluded) more than chance would produce?

Two null models are run:
  A (naive)      -- per dataset, draw n_j significant circRNAs uniformly at
                    random from that dataset's tested universe.
  B (stratified, main result) -- circRNAs are stratified by breadth of
                    testing (how many datasets test them at all) into 5
                    quantile strata. Within each dataset, the real
                    significant list's distribution across strata is
                    preserved when drawing the null replacement set, so
                    "easy to detect everywhere" circRNAs remain
                    proportionally more likely to be drawn in the null too.

circRNA identity across datasets = exact string match on circ_id
("chr:start|end", pipe-delimited, no strand field -- this pipeline's native
format), consistent with web_ui.py's _load_cross_dataset_data(), which does
the same exact-match join for the live /cross_dataset recurrence page.

Usage:
  python run_permutation.py --out-dir results/permutation \
      [--exclude GSE192410] [--n-perm 10000] [--seed 20260930]
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

BASE_DIR = Path(os.path.expanduser("~")) / "circRNA_agent"

ALL_19 = [
    "GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998",
    "GSE136569", "GSE143797", "GSE148036", "GSE171011", "GSE192410",
    "GSE192849", "GSE221107", "GSE229705", "GSE248612", "GSE58135",
    "GSE77509", "GSE97239", "PRJNA553289", "SRP156355",
]  # GSE323364 (cell line, not tumor/normal) excluded from the main analysis


# ── Data loading ────────────────────────────────────────────────────────────

def load_tested_and_sig(gse_list: list[str]) -> tuple[dict, dict, dict]:
    """Returns (tested[gse] -> set(circ_id), sig[gse] -> set(circ_id), n_sig[gse] -> int)."""
    tested, sig, n_sig_known = {}, {}, {}
    for gse in gse_list:
        cfg_path = BASE_DIR / "config" / "projects" / f"{gse}.yaml"
        cfg = yaml.safe_load(cfg_path.read_text())
        rdir = Path(cfg["results_dir"])
        de_cfg = cfg.get("de", {})
        fdr_thr = float(de_cfg.get("fdr_cutoff", 0.05))
        lfc_thr = float(de_cfg.get("log2fc_cutoff", 1.0))
        sig_by = de_cfg.get("de_sig_by", "pvalue")
        p_col = "pvalue" if sig_by == "pvalue" else "padj"

        df = pd.read_csv(rdir / "de" / "de_results.tsv", sep="\t", low_memory=False)
        df["circ_id"] = df["circ_id"].astype(str)
        tested[gse] = set(df["circ_id"])
        mask = (df[p_col] < fdr_thr) & (df["log2FC"].abs() > lfc_thr)
        sig[gse] = set(df.loc[mask, "circ_id"])
        n_sig_known[gse] = int(mask.sum())
    return tested, sig, n_sig_known


# ── Step 1: matrices + validation ───────────────────────────────────────────

def build_matrices(gse_list: list[str], tested: dict, sig: dict):
    all_ids = sorted(set().union(*tested.values()))
    id_index = {cid: i for i, cid in enumerate(all_ids)}
    n_circ, n_ds = len(all_ids), len(gse_list)

    T = np.zeros((n_circ, n_ds), dtype=bool)
    S = np.zeros((n_circ, n_ds), dtype=bool)
    for j, gse in enumerate(gse_list):
        for cid in tested[gse]:
            T[id_index[cid], j] = True
        for cid in sig[gse]:
            S[id_index[cid], j] = True
    return all_ids, T, S


def validate_matrices(gse_list, T, S, n_sig_known) -> bool:
    ok = True
    col_sums_S = S.sum(axis=0)
    for j, gse in enumerate(gse_list):
        if col_sums_S[j] != n_sig_known[gse]:
            print(f"[VALIDATION FAIL] {gse}: S colsum={col_sums_S[j]} != known n_sig={n_sig_known[gse]}")
            ok = False
    if not np.all(S <= T):
        n_bad = int(np.sum(S & ~T))
        print(f"[VALIDATION FAIL] S is not a subset of T: {n_bad} entries significant but not tested")
        ok = False
    if ok:
        print("[VALIDATION OK] S column sums match known significant counts; S subset-of T holds.")
    return ok


# ── Step 2: observed recurrence ─────────────────────────────────────────────

def observed_recurrence(S: np.ndarray) -> np.ndarray:
    return S.sum(axis=1)


def recurrence_summary(rec: np.ndarray) -> dict:
    return {
        "n>=2": int((rec >= 2).sum()),
        "n>=3": int((rec >= 3).sum()),
        "n>=4": int((rec >= 4).sum()),
        "n>=5": int((rec >= 5).sum()),
        "max_recurrence": int(rec.max()) if len(rec) else 0,
    }


# ── Step 3: null models ─────────────────────────────────────────────────────

def null_naive(T: np.ndarray, n_sig_per_ds: np.ndarray, n_perm: int, rng) -> np.ndarray:
    """Version A: per dataset, draw n_j significant circRNAs uniformly from T[:,j]."""
    n_circ, n_ds = T.shape
    results = np.zeros((n_perm, 5), dtype=np.int64)  # [n>=2,n>=3,n>=4,n>=5,max]
    tested_idx = [np.where(T[:, j])[0] for j in range(n_ds)]
    for p in range(n_perm):
        rec = np.zeros(n_circ, dtype=np.int32)
        for j in range(n_ds):
            n_j = n_sig_per_ds[j]
            if n_j == 0 or len(tested_idx[j]) == 0:
                continue
            draw = rng.choice(tested_idx[j], size=min(n_j, len(tested_idx[j])), replace=False)
            rec[draw] += 1
        results[p] = [(rec >= 2).sum(), (rec >= 3).sum(), (rec >= 4).sum(),
                      (rec >= 5).sum(), rec.max()]
    return results


def assign_strata(T: np.ndarray, n_strata: int = 5) -> np.ndarray:
    """Stratify circRNAs by breadth of testing (rowSums(T)) into n_strata quantile bins."""
    breadth = T.sum(axis=1)
    # Quantile-based bin edges on the *tested* circRNAs (breadth >= 1)
    nz = breadth[breadth > 0]
    edges = np.unique(np.quantile(nz, np.linspace(0, 1, n_strata + 1)))
    strata = np.digitize(breadth, edges[1:-1], right=True)
    strata[breadth == 0] = -1  # never tested anywhere; shouldn't occur given construction
    return strata


def null_stratified(T: np.ndarray, S: np.ndarray, strata: np.ndarray,
                     n_perm: int, rng) -> np.ndarray:
    """Version B: preserve the real significant list's per-stratum proportion
    within each dataset when drawing the null replacement set."""
    n_circ, n_ds = T.shape
    n_strata = int(strata.max()) + 1
    results = np.zeros((n_perm, 5), dtype=np.int64)

    # Precompute, per dataset, per stratum: tested-pool indices, and the
    # observed number of real significant circRNAs in that stratum (to set
    # the draw size for the null).
    tested_pool = {}
    draw_size = {}
    for j in range(n_ds):
        for k in range(n_strata):
            pool = np.where(T[:, j] & (strata == k))[0]
            tested_pool[(j, k)] = pool
            draw_size[(j, k)] = int((S[:, j] & (strata == k)).sum())

    for p in range(n_perm):
        rec = np.zeros(n_circ, dtype=np.int32)
        for j in range(n_ds):
            for k in range(n_strata):
                pool = tested_pool[(j, k)]
                n_draw = draw_size[(j, k)]
                if n_draw == 0 or len(pool) == 0:
                    continue
                draw = rng.choice(pool, size=min(n_draw, len(pool)), replace=False)
                rec[draw] += 1
        results[p] = [(rec >= 2).sum(), (rec >= 3).sum(), (rec >= 4).sum(),
                      (rec >= 5).sum(), rec.max()]
    return results


# ── Step 4: stats ────────────────────────────────────────────────────────────

def summarize_null(null_arr: np.ndarray, observed: dict, label: str) -> list[dict]:
    cols = ["n>=2", "n>=3", "n>=4", "n>=5", "max_recurrence"]
    rows = []
    n_perm = null_arr.shape[0]
    for i, col in enumerate(cols):
        null_col = null_arr[:, i]
        obs = observed[col]
        n_ge = int((null_col >= obs).sum())
        p_emp = (1 + n_ge) / (1 + n_perm)
        null_median = float(np.median(null_col))
        rows.append({
            "null_version": label,
            "metric": col,
            "observed": obs,
            "null_median": null_median,
            "null_mean": float(np.mean(null_col)),
            "null_p2.5": float(np.percentile(null_col, 2.5)),
            "null_p97.5": float(np.percentile(null_col, 97.5)),
            "empirical_p": f"<{1/(1+n_perm):.1e}" if n_ge == 0 else f"{p_emp:.4g}",
            "effect_size_vs_median": (obs / null_median) if null_median > 0 else float("inf"),
        })
    return rows


# ── Main ─────────────────────────────────────────────────────────────────────

def run_analysis(gse_list: list[str], label: str, n_perm: int, seed: int, out_dir: Path):
    print(f"\n{'='*70}\n{label}  (n_datasets={len(gse_list)})\n{'='*70}")
    tested, sig, n_sig_known = load_tested_and_sig(gse_list)
    all_ids, T, S = build_matrices(gse_list, tested, sig)
    ok = validate_matrices(gse_list, T, S, n_sig_known)
    if not ok:
        raise SystemExit(f"[{label}] Validation failed -- stopping, see messages above.")

    rec = observed_recurrence(S)
    obs_summary = recurrence_summary(rec)
    print(f"[{label}] Observed: {obs_summary}")

    n_sig_per_ds = S.sum(axis=0)
    rng_a = np.random.default_rng(seed)
    rng_b = np.random.default_rng(seed)

    null_a = null_naive(T, n_sig_per_ds, n_perm, rng_a)
    strata = assign_strata(T, n_strata=5)
    null_b = null_stratified(T, S, strata, n_perm, rng_b)

    rows_a = summarize_null(null_a, obs_summary, "A_naive")
    rows_b = summarize_null(null_b, obs_summary, "B_stratified")
    for r in rows_a + rows_b:
        r["dataset_set"] = label
        r["n_datasets"] = len(gse_list)
        r["n_perm"] = n_perm
        r["seed"] = seed

    return rows_a + rows_b, null_b, obs_summary, rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(BASE_DIR / "results" / "permutation"))
    ap.add_argument("--exclude", default="GSE192410",
                     help="Dataset to drop for the sensitivity analysis")
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260930)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []

    # Main analysis: all 19 (GSE323364 cell-line already excluded from ALL_19)
    rows_main, null_b_main, obs_main, rec_main = run_analysis(
        ALL_19, "main_19_datasets", args.n_perm, args.seed, out_dir)
    all_rows += rows_main

    # Sensitivity: drop the specified dataset
    sens_list = [g for g in ALL_19 if g != args.exclude]
    assert len(sens_list) == len(ALL_19) - 1, f"{args.exclude} not found in ALL_19"
    rows_sens, null_b_sens, obs_sens, rec_sens = run_analysis(
        sens_list, f"sensitivity_excl_{args.exclude}", args.n_perm, args.seed, out_dir)
    all_rows += rows_sens

    df = pd.DataFrame(all_rows)
    tsv_path = out_dir / "permutation_summary.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"\n[OUT] {tsv_path}")

    # ── Plot (version B, main analysis) ──────────────────────────────────────
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    metrics = ["n>=2", "n>=3", "n>=4", "n>=5"]
    col_idx = {"n>=2": 0, "n>=3": 1, "n>=4": 2, "n>=5": 3}
    for ax, metric in zip(axes.flat, metrics):
        null_vals = null_b_main[:, col_idx[metric]]
        obs_val = obs_main[metric]
        row = df[(df["dataset_set"] == "main_19_datasets") &
                 (df["null_version"] == "B_stratified") &
                 (df["metric"] == metric)].iloc[0]
        ax.hist(null_vals, bins=40, color="#94a3b8", edgecolor="white")
        ax.axvline(obs_val, color="#dc2626", linewidth=2)
        ax.set_title(f"{metric}: observed={obs_val}, p={row['empirical_p']}", fontsize=11)
        ax.set_xlabel("# circRNAs under null (stratified)")
        ax.set_ylabel("permutation count")
    fig.suptitle(f"Cross-dataset recurrence vs. stratified null (n_perm={args.n_perm}, "
                 f"n_datasets=19, seed={args.seed})", fontsize=12)
    fig.tight_layout()
    png_path = out_dir / "permutation_null.png"
    fig.savefig(png_path, dpi=150)
    print(f"[OUT] {png_path}")

    # ── Methods paragraph ────────────────────────────────────────────────────
    main_b = df[(df["dataset_set"] == "main_19_datasets") & (df["null_version"] == "B_stratified")]
    p_ge2 = main_b[main_b["metric"] == "n>=2"]["empirical_p"].iloc[0]
    p_ge5 = main_b[main_b["metric"] == "n>=5"]["empirical_p"].iloc[0]
    obs_ge2 = main_b[main_b["metric"] == "n>=2"]["observed"].iloc[0]
    obs_ge5 = main_b[main_b["metric"] == "n>=5"]["observed"].iloc[0]
    med_ge2 = main_b[main_b["metric"] == "n>=2"]["null_median"].iloc[0]
    med_ge5 = main_b[main_b["metric"] == "n>=5"]["null_median"].iloc[0]

    zh = f"""## Methods：跨資料集重現性之置換檢定

為評估觀察到的跨資料集 circRNA 重現數（{obs_main['n>=2']} / {obs_main['n>=3']} / {obs_main['n>=4']} / {obs_main['n>=5']} 個 circRNA 分別於 ≥2 / ≥3 / ≥4 / ≥5 個資料集同時達顯著，共 19 個 tumor-vs-normal 資料集，排除非腫瘤對比之細胞株資料集 GSE323364）是否顯著高於隨機期望，對每個資料集分別以其實際顯著 circRNA 數 $n_j$ 為抽樣數，從該資料集之「可測 circRNA 全集」（filterByExpr 篩選後、DE 檢定前之完整清單，而非僅顯著清單）中重新隨機指派顯著標籤，建立虛無分布。circRNA 身分跨資料集以基因體座標（`chr:start|end`）完全字串相等比對，與既有跨資料集重現性分析（`/cross_dataset` 頁面）所用規則一致。

為避免簡單隨機抽樣低估虛無重現數（因偵測頻率高、表現量穩定的 circRNA 本來就較易在多個資料集中被判定顯著），採用分層置換法：依各 circRNA 被納入檢定之資料集數（testing breadth）分為 5 個分位層，於各資料集內依其真實顯著清單在各層之分布比例，自對應層之可測集合中抽樣，使虛無分布保留「廣泛可測之 circRNA 較易中選」此一結構性偏誤。另以簡單隨機抽樣（不分層）作為上界對照。兩版本均執行 10,000 次置換（`numpy.random.default_rng(20260930)` 固定亂數種子），經驗 p 值計算為 $p = (1+\\#\\{{\\text{{null}} \\geq \\text{{observed}}\\}})/(1+10000)$。

結果顯示，在分層虛無模型下，≥2 個資料集重現的虛無分布中位數為 {med_ge2:.1f}（觀察值 {obs_ge2}，p={p_ge2}），≥5 個資料集重現的虛無分布中位數為 {med_ge5:.1f}（觀察值 {obs_ge5}，p={p_ge5}）。為評估 GSE192410（因 DCC 偵測方法不對稱已修正，詳見方法學章節）對結果之影響，另排除該資料集重複整組分析作為敏感度分析，結果列於補充表。"""

    en = f"""## Methods: Permutation Test for Cross-Dataset Recurrence

To assess whether the observed cross-dataset recurrence of significant circRNAs ({obs_main['n>=2']} / {obs_main['n>=3']} / {obs_main['n>=4']} / {obs_main['n>=5']} circRNAs reaching significance in ≥2 / ≥3 / ≥4 / ≥5 of 19 tumor-vs-normal datasets, with the non-tumor cell-line dataset GSE323364 excluded) exceeds chance expectation, we built a null distribution by re-assigning, independently per dataset, $n_j$ "significant" labels (matching the dataset's true significant count) uniformly at random among that dataset's tested circRNA universe (the full post-filterByExpr, pre-significance-filtering set, not merely the significant list). Cross-dataset circRNA identity was resolved by exact string match on genomic coordinates (`chr:start|end`), consistent with the matching rule used by the existing cross-dataset recurrence page (`/cross_dataset`).

Because naive uniform resampling under-estimates null recurrence (circRNAs tested across many datasets, typically those with stable, detectable expression, are intrinsically more likely to be called significant anywhere), we additionally ran a stratified permutation: circRNAs were binned into 5 quantile strata by testing breadth (number of datasets in which they were tested at all), and within each dataset, the real significant list's per-stratum proportions were preserved when drawing the null replacement set — so the null retains the structural bias that broadly-tested circRNAs are easier to call significant. A naive (unstratified) permutation was also run as an upper-bound reference. Both versions used 10,000 permutations (`numpy.random.default_rng(20260930)`), with empirical p-values computed as $p = (1+\\#\\{{\\text{{null}} \\geq \\text{{observed}}\\}})/(1+10000)$.

Under the stratified null, the median null count for ≥2-dataset recurrence was {med_ge2:.1f} (observed {obs_ge2}, p={p_ge2}), and for ≥5-dataset recurrence {med_ge5:.1f} (observed {obs_ge5}, p={p_ge5}). To assess the influence of GSE192410 (whose DCC-detection asymmetry was corrected; see Methods), the full analysis was repeated excluding this dataset as a sensitivity check, reported in Supplementary Table."""

    md_path = out_dir / "methods_paragraph.md"
    md_path.write_text(zh + "\n\n---\n\n" + en, encoding="utf-8")
    print(f"[OUT] {md_path}")

    print("\n[DONE] All outputs written to", out_dir)


if __name__ == "__main__":
    main()
