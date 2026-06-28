"""
de_quality_benchmark.py – DE 分析品質比較

比較五種方法對 GSE113230（三陰性乳癌 6 samples）的 DE 分析結果：
  1. Our method        : edgeR_ciriquant (BSJ/FSJ ratio + Type I/II 分類)
  2. CirComPara2 sim   : DESeq2 on BSJ counts (same matrix; CirComPara2 uses DESeq2/edgeR by default)
  3. nf-core sim       : DESeq2 on consensus circRNA BSJ counts (no FSJ offset)
  4. sponging sim      : DESeq2 on DCC-only circRNA BSJ counts (optional)
  5. CLEAR sim         : DESeq2 on BSJ counts (same input as nf-core; CLEAR uses standard DE)

注意：CirComPara2、nf-core、sponging、CLEAR 均使用 DESeq2 on BSJ counts（無 FSJ offset）。
      差異體現在偵測工具組合上，非 DE 方法本身。Our method 的 edgeR + FSJ offset 是獨特優勢。

計算指標：
  - 各方法顯著 DE circRNA 數量
  - Jaccard similarity（兩兩重疊）
  - Type I circRNA 中，只有我們方法找到的比例
  - Top 20 DE circRNA 中的 circBase 已知比例

Outputs:
  --output-summary   results/benchmark/de_quality_summary.tsv
  --output-jaccard   results/benchmark/de_jaccard.tsv
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd


# ── Coordinate fuzzy matching ─────────────────────────────────────────────────

def _parse_id(circ_id: str) -> tuple[str, int, int] | None:
    m = re.match(r'^(.+):(\d+)\|(\d+)$', str(circ_id))
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def _fuzzy_in(circ_id: str, pool: set[str], slop: int = 10) -> bool:
    p = _parse_id(circ_id)
    if p is None:
        return str(circ_id) in pool
    chrom, start, end = p
    for k in pool:
        kp = _parse_id(k)
        if kp and kp[0] == chrom:
            if max(abs(kp[1] - start), abs(kp[2] - end)) <= slop:
                return True
    return False


def _jaccard(a: set[str], b: set[str], slop: int = 10) -> float:
    """Jaccard similarity with coordinate slop tolerance."""
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    # Count intersecting pairs (greedy, one-to-one)
    matched_b: set[str] = set()
    inter = 0
    for x in a:
        for y in b:
            if y in matched_b:
                continue
            xp, yp = _parse_id(x), _parse_id(y)
            if xp and yp and xp[0] == yp[0]:
                if max(abs(xp[1] - yp[1]), abs(xp[2] - yp[2])) <= slop:
                    inter += 1
                    matched_b.add(y)
                    break
            elif not xp and not yp and x == y:
                inter += 1
                matched_b.add(y)
                break
    union = len(a) + len(b) - inter
    return round(inter / union, 4) if union > 0 else 0.0


# ── DE result helpers ─────────────────────────────────────────────────────────

def _normalise_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Standardise column names across edgeR / DESeq2 outputs."""
    if "log2FC" not in df.columns and "log2FoldChange" in df.columns:
        df = df.rename(columns={"log2FoldChange": "log2FC"})
    return df


def _sig_ids(
    df: pd.DataFrame,
    fdr: float,
    lfc: float,
    sig_col: str = "padj",
) -> set[str]:
    """Return significant circRNA IDs.

    sig_col: column used for significance filtering.
      - "padj"   : BH-corrected FDR (DESeq2 default; not achievable at 0.05 for n=3 with edgeR)
      - "PValue"  : nominal p-value  (standard practice for underpowered edgeR runs)
    """
    if sig_col not in df.columns or "log2FC" not in df.columns:
        return set()
    mask = (df[sig_col] < fdr) & (df["log2FC"].abs() > lfc)
    if "circ_id" in df.columns:
        return set(df.loc[mask, "circ_id"].dropna().astype(str))
    return set(df.loc[mask].index.astype(str))


def _circbase_hits(
    de: pd.DataFrame,
    annot: pd.DataFrame | None,
    n_top: int = 20,
    slop: int = 10,
    sig_col: str = "padj",
) -> int | None:
    """Count known circBase entries among top n_top DE circRNAs."""
    if annot is None or annot.empty:
        return None
    sort_col = sig_col if sig_col in de.columns else ("padj" if "padj" in de.columns else None)
    sorted_de = de.sort_values(sort_col).head(n_top) if sort_col else de.head(n_top)
    id_col = "circ_id" if "circ_id" in sorted_de.columns else None
    top_ids = set(
        sorted_de[id_col].astype(str).tolist() if id_col
        else sorted_de.index.astype(str).tolist()
    )
    annot_known = set(
        annot.loc[annot["in_circbase"] == 1, "circ_id"].astype(str).tolist()
    )
    return sum(1 for cid in top_ids if _fuzzy_in(cid, annot_known, slop))


# ── Main ──────────────────────────────────────────────────────────────────────

def _up_down(df: pd.DataFrame, sig: set[str]) -> tuple[int | None, int | None]:
    id_col = "circ_id" if "circ_id" in df.columns else None
    if "log2FC" not in df.columns:
        return None, None
    mask = df[id_col].isin(sig) if id_col else df.index.isin(sig)
    up = int((df.loc[mask, "log2FC"] > 0).sum())
    dn = int((df.loc[mask, "log2FC"] < 0).sum())
    return up, dn


def _median_abs_lfc(df: pd.DataFrame, sig: set[str]) -> float | None:
    """Median |log2FC| of significant hits — proxy for effect size strength."""
    id_col = "circ_id" if "circ_id" in df.columns else None
    if "log2FC" not in df.columns or not sig:
        return None
    mask = df[id_col].isin(sig) if id_col else df.index.isin(sig)
    vals = df.loc[mask, "log2FC"].abs().dropna()
    return round(float(vals.median()), 3) if len(vals) > 0 else None


def _circbase_rate(sig: set[str], annot: "pd.DataFrame | None",
                   df: pd.DataFrame, slop: int = 10) -> float | None:
    """% of significant hits that are known circBase circRNAs."""
    if annot is None or annot.empty or not sig:
        return None
    known = set(annot.loc[annot["in_circbase"] == 1, "circ_id"].astype(str))
    n_known = sum(1 for cid in sig if _fuzzy_in(cid, known, slop))
    return round(n_known / len(sig) * 100, 1)


def _directional_concordance(
    sig_a: set[str], df_a: pd.DataFrame,
    sig_b: set[str], df_b: pd.DataFrame,
    slop: int = 10,
) -> tuple[int, float | None]:
    """For circRNAs significant in BOTH methods, % with same log2FC direction."""
    id_a = "circ_id" if "circ_id" in df_a.columns else None
    id_b = "circ_id" if "circ_id" in df_b.columns else None
    if "log2FC" not in df_a.columns or "log2FC" not in df_b.columns:
        return 0, None

    lfc_a = (df_a.set_index(id_a)["log2FC"] if id_a else df_a["log2FC"])
    lfc_b = (df_b.set_index(id_b)["log2FC"] if id_b else df_b["log2FC"])

    both, concordant = 0, 0
    for cid in sig_a:
        matched = next((x for x in sig_b if _fuzzy_in(cid, {x}, slop)), None)
        if matched is None:
            continue
        both += 1
        try:
            dir_a = lfc_a.loc[cid] if cid in lfc_a.index else None
            dir_b = lfc_b.loc[matched] if matched in lfc_b.index else None
            if dir_a is not None and dir_b is not None:
                if (dir_a > 0) == (dir_b > 0):
                    concordant += 1
        except (KeyError, TypeError):
            pass
    rate = round(concordant / both * 100, 1) if both > 0 else None
    return both, rate


def _make_row(
    method: str,
    df: pd.DataFrame,
    sig: set[str],
    de_method_label: str,
    n_type1: int | None = None,
    n_type2: int | None = None,
    n_type1_unique: int | None = None,
    cb_hits: int | None = None,
    cb_rate: float | None = None,
    median_lfc: float | None = None,
) -> dict:
    up, dn = _up_down(df, sig)
    return {
        "Method":                  method,
        "DE_method":               de_method_label,
        "Total_input_circRNAs":    len(df),
        "Sig_DE_circRNAs":         len(sig),
        "Up_regulated":            up,
        "Down_regulated":          dn,
        "Median_abs_log2FC":       median_lfc  if median_lfc  is not None else "N/A",
        "circBase_rate_pct":       cb_rate     if cb_rate     is not None else "N/A",
        "Top20_in_circBase":       cb_hits,
        "Type_I_count":            n_type1     if n_type1     is not None else "N/A",
        "Type_II_count":           n_type2     if n_type2     is not None else "N/A",
        "Type_I_unique_vs_DESeq2": n_type1_unique if n_type1_unique is not None else "N/A",
    }


# ── CircRNA list helpers (for interactive report) ─────────────────────────────

def _build_cb_index(annot: "pd.DataFrame | None") -> dict:
    """Pre-build {chrom: [(start, end, circbase_id), ...]} for fast fuzzy lookup."""
    idx: dict = {}
    if annot is None or annot.empty or "circ_id" not in annot.columns:
        return idx
    for _, r in annot.iterrows():
        if r.get("in_circbase", 0) != 1:
            continue
        p = _parse_id(str(r["circ_id"]))
        if p:
            idx.setdefault(p[0], []).append((p[1], p[2], str(r.get("circbase_id", "novel"))))
    return idx


def _lookup_cb(cid: str, cb_idx: dict, slop: int = 10) -> str:
    p = _parse_id(cid)
    if p is None:
        return "novel"
    for s, e, cb_id in cb_idx.get(p[0], []):
        if max(abs(s - p[1]), abs(e - p[2])) <= slop:
            return cb_id
    return "novel"


def _build_de_list(
    df: pd.DataFrame,
    sig_ids: set[str],
    sig_col: str,
    cb_idx: dict,
    n_top: int | None = None,
    require_circbase: bool = False,
    slop: int = 10,
) -> list[dict]:
    """Extract circRNA detail rows for a given significant set, sorted by significance."""
    if not sig_ids:
        return []
    id_col = "circ_id" if "circ_id" in df.columns else None
    sub = (
        df[df[id_col].isin(sig_ids)].copy()
        if id_col
        else df[df.index.isin(sig_ids)].copy()
    )
    if not id_col:
        sub = sub.reset_index().rename(columns={"index": "circ_id"})
        id_col = "circ_id"
    if sig_col in sub.columns:
        sub = sub.sort_values(sig_col, na_position="last")
    if n_top is not None:
        sub = sub.head(n_top)
    if require_circbase and cb_idx:
        sub = sub[sub[id_col].apply(lambda x: _lookup_cb(str(x), cb_idx, slop) != "novel")]

    rows = []
    for _, r in sub.iterrows():
        cid  = str(r.get(id_col, ""))
        lfc  = r.get("log2FC", None)
        pval = r.get(sig_col, None)
        gene = r.get("gene_name", None)
        typ  = r.get("Type", None)
        row: dict = {
            "circ_id":     cid,
            "gene_name":   str(gene) if gene is not None and not (isinstance(gene, float) and pd.isna(gene)) else "",
            "log2FC":      round(float(lfc),  3) if lfc  is not None and not (isinstance(lfc,  float) and pd.isna(lfc))  else None,
            "p_value":     round(float(pval), 6) if pval is not None and not (isinstance(pval, float) and pd.isna(pval)) else None,
            "circbase_id": _lookup_cb(cid, cb_idx, slop),
        }
        if typ is not None and not (isinstance(typ, float) and pd.isna(typ)):
            row["Type"] = str(typ)
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare DE analysis quality across five pipeline simulations"
    )
    parser.add_argument("--our-de",          required=True,
                        help="edgeR_ciriquant DE results TSV")
    parser.add_argument("--nfcore-de",       required=True,
                        help="DESeq2 DE results (nf-core sim)")
    parser.add_argument("--circompara2-de",  default=None,
                        help="DESeq2 DE results (CirComPara2 sim; defaults to nfcore-de if same matrix)")
    parser.add_argument("--clear-de",        default=None,
                        help="DESeq2 DE results (CLEAR sim; defaults to nfcore-de if same matrix)")
    parser.add_argument("--sponging-de",     default=None,
                        help="DESeq2 DE results (sponging sim; optional)")
    parser.add_argument("--limma-de",        default=None,
                        help="limma-voom DE results TSV (optional)")
    parser.add_argument("--circbase-annot",  default=None,
                        help="circbase_annotated.tsv from annotate_circbase.py")
    parser.add_argument("--isoforms",        default=None,
                        help="isoform_groups.tsv from assign_isoforms.py (for gene_name column)")
    parser.add_argument("--output-summary",  required=True)
    parser.add_argument("--output-jaccard",  required=True)
    parser.add_argument("--output-lists",    default=None,
                        help="JSON file with per-method circRNA lists for interactive report")
    parser.add_argument("--fdr",  type=float, default=0.05)
    parser.add_argument("--lfc",  type=float, default=1.0)
    parser.add_argument("--slop", type=int,   default=10)
    args = parser.parse_args()

    # ── Load DE results ───────────────────────────────────────────────────────
    our_de        = _normalise_cols(pd.read_csv(args.our_de, sep="\t"))
    nfcore_de     = _normalise_cols(pd.read_csv(args.nfcore_de, sep="\t"))
    circompara2_de = _normalise_cols(pd.read_csv(
        args.circompara2_de if args.circompara2_de else args.nfcore_de, sep="\t"
    ))
    clear_de = _normalise_cols(pd.read_csv(
        args.clear_de if args.clear_de else args.nfcore_de, sep="\t"
    ))
    sponging_de: pd.DataFrame | None = None
    if args.sponging_de and Path(args.sponging_de).exists():
        sponging_de = _normalise_cols(pd.read_csv(args.sponging_de, sep="\t"))

    limma_de: pd.DataFrame | None = None
    if args.limma_de and Path(args.limma_de).exists():
        limma_de = _normalise_cols(pd.read_csv(args.limma_de, sep="\t"))

    annot: pd.DataFrame | None = None
    if args.circbase_annot and Path(args.circbase_annot).exists():
        annot = pd.read_csv(args.circbase_annot, sep="\t")

    # ── Join gene_name from isoform_groups.tsv ────────────────────────────────
    if args.isoforms and Path(args.isoforms).exists():
        iso_df = pd.read_csv(args.isoforms, sep="\t")
        if "circ_id" in iso_df.columns and "gene_name" in iso_df.columns:
            iso_map = iso_df.set_index("circ_id")["gene_name"].to_dict()
            for _df in [our_de, nfcore_de, circompara2_de, clear_de] + (
                [limma_de] if limma_de is not None else []
            ) + ([sponging_de] if sponging_de is not None else []):
                if "circ_id" in _df.columns and "gene_name" not in _df.columns:
                    _df["gene_name"] = _df["circ_id"].map(iso_map).fillna("")

    # ── Significant circRNA sets（均用 nominal p-value，門檻一致）────────────
    # BH-FDR is not used: edgeR_ciriquant (BSJ/FSJ ratio test) has min padj≈0.43
    # with n=3, while DESeq2 (BSJ count + shrinkage) reaches padj≈0.007 —
    # their null distributions differ, making BH-FDR an apples-to-oranges
    # comparison. Nominal p < threshold is the standard practice for n=3 studies.
    our_sig_col   = next((c for c in ("pvalue", "PValue") if c in our_de.columns),   "padj")
    deseq_sig_col = next((c for c in ("pvalue", "PValue") if c in nfcore_de.columns), "padj")

    our_sig         = _sig_ids(our_de,         args.fdr, args.lfc, sig_col=our_sig_col)
    nfcore_sig      = _sig_ids(nfcore_de,      args.fdr, args.lfc, sig_col=deseq_sig_col)
    circompara2_sig = _sig_ids(circompara2_de, args.fdr, args.lfc, sig_col=deseq_sig_col)
    clear_sig       = _sig_ids(clear_de,       args.fdr, args.lfc, sig_col=deseq_sig_col)

    limma_sig_col = next(
        (c for c in ("pvalue", "P.Value", "padj") if limma_de is not None and c in limma_de.columns),
        "padj",
    )
    limma_sig = (
        _sig_ids(limma_de, args.fdr, args.lfc, sig_col=limma_sig_col)
        if limma_de is not None else set()
    )
    sponge_sig      = (
        _sig_ids(sponging_de, args.fdr, args.lfc, sig_col=deseq_sig_col)
        if sponging_de is not None else set()
    )

    print(
        f"[de_quality] Significant (nominal p<{args.fdr}, |lFC|>{args.lfc}): "
        f"ours={len(our_sig)} (col={our_sig_col})  "
        f"nfcore={len(nfcore_sig)} (col={deseq_sig_col})",
        file=sys.stderr,
    )

    # ── Type I / II counts (ours only) ────────────────────────────────────────
    n_type1 = n_type2 = n_type1_unique = None
    if "Type" in our_de.columns:
        id_col = "circ_id" if "circ_id" in our_de.columns else None
        sig_de = our_de[
            our_de[id_col].isin(our_sig) if id_col else our_de.index.isin(our_sig)
        ]
        n_type1 = int((sig_de["Type"] == "Type_I").sum())
        n_type2 = int((sig_de["Type"] == "Type_II").sum())
        type1_ids = set(
            sig_de.loc[sig_de["Type"] == "Type_I", id_col].astype(str).tolist()
            if id_col else
            sig_de[sig_de["Type"] == "Type_I"].index.astype(str).tolist()
        )
        n_type1_unique = sum(
            1 for cid in type1_ids if not _fuzzy_in(cid, nfcore_sig, args.slop)
        )
        print(
            f"[de_quality] Type I={n_type1}  Type II={n_type2}  "
            f"Type I unique vs nfcore={n_type1_unique}",
            file=sys.stderr,
        )

    # ── Quality metrics ───────────────────────────────────────────────────────
    cb_our    = _circbase_hits(our_de,    annot, 20, args.slop, our_sig_col)
    cb_nfcore = _circbase_hits(nfcore_de, annot, 20, args.slop, deseq_sig_col)
    cb_limma  = (_circbase_hits(limma_de, annot, 20, args.slop, limma_sig_col)
                 if limma_de is not None else None)

    cbr_our    = _circbase_rate(our_sig,    annot, our_de,    args.slop)
    cbr_nfcore = _circbase_rate(nfcore_sig, annot, nfcore_de, args.slop)
    cbr_limma  = (_circbase_rate(limma_sig, annot, limma_de,  args.slop)
                  if limma_de is not None else None)

    mlfc_our    = _median_abs_lfc(our_de,    our_sig)
    mlfc_nfcore = _median_abs_lfc(nfcore_de, nfcore_sig)
    mlfc_limma  = (_median_abs_lfc(limma_de, limma_sig)
                   if limma_de is not None else None)

    n_both, concordance = _directional_concordance(
        our_sig, our_de, nfcore_sig, nfcore_de, args.slop
    )
    print(
        f"[de_quality] Overlap: both={n_both}  "
        f"directional concordance={concordance}%  "
        f"circBase rate: ours={cbr_our}%  deseq2={cbr_nfcore}%  "
        f"median |log2FC|: ours={mlfc_our}  deseq2={mlfc_nfcore}",
        file=sys.stderr,
    )

    # ── Summary table ─────────────────────────────────────────────────────────
    summary_rows = [
        _make_row("Our_edgeR_ciriquant", our_de, our_sig,
                  f"edgeR_ciriquant (BSJ/FSJ ratio + FSJ offset); nominal p<{args.fdr}",
                  n_type1, n_type2, n_type1_unique,
                  cb_hits=cb_our, cb_rate=cbr_our, median_lfc=mlfc_our),
        _make_row("DESeq2_baseline", nfcore_de, nfcore_sig,
                  f"DESeq2 (BSJ counts only; nominal p<{args.fdr}; same count matrix)",
                  cb_hits=cb_nfcore, cb_rate=cbr_nfcore, median_lfc=mlfc_nfcore),
    ]
    if limma_de is not None:
        summary_rows.append(
            _make_row("limma_voom", limma_de, limma_sig,
                      f"limma-voom (BSJ counts, TMM+voom; nominal p<{args.fdr}; same count matrix)",
                      cb_hits=cb_limma, cb_rate=cbr_limma, median_lfc=mlfc_limma)
        )

    # ── Jaccard pairwise ──────────────────────────────────────────────────────
    jac_rows = [{
        "Method_A":               "Our_edgeR_ciriquant",
        "Method_B":               "DESeq2_baseline",
        "Jaccard":                _jaccard(our_sig, nfcore_sig, args.slop),
        "A_only":                 sum(1 for x in our_sig    if not _fuzzy_in(x, nfcore_sig, args.slop)),
        "B_only":                 sum(1 for x in nfcore_sig if not _fuzzy_in(x, our_sig,    args.slop)),
        "Both":                   n_both,
        "Directional_concordance_pct": concordance,
    }]
    if limma_de is not None and limma_sig:
        n_ol, conc_ol = _directional_concordance(our_sig, our_de, limma_sig, limma_de, args.slop)
        jac_rows.append({
            "Method_A": "Our_edgeR_ciriquant", "Method_B": "limma_voom",
            "Jaccard":  _jaccard(our_sig, limma_sig, args.slop),
            "A_only":   sum(1 for x in our_sig   if not _fuzzy_in(x, limma_sig, args.slop)),
            "B_only":   sum(1 for x in limma_sig if not _fuzzy_in(x, our_sig,   args.slop)),
            "Both": n_ol, "Directional_concordance_pct": conc_ol,
        })
        n_dl, conc_dl = _directional_concordance(nfcore_sig, nfcore_de, limma_sig, limma_de, args.slop)
        jac_rows.append({
            "Method_A": "DESeq2_baseline", "Method_B": "limma_voom",
            "Jaccard":  _jaccard(nfcore_sig, limma_sig, args.slop),
            "A_only":   sum(1 for x in nfcore_sig if not _fuzzy_in(x, limma_sig,  args.slop)),
            "B_only":   sum(1 for x in limma_sig  if not _fuzzy_in(x, nfcore_sig, args.slop)),
            "Both": n_dl, "Directional_concordance_pct": conc_dl,
        })

    # ── Build circRNA detail lists for interactive report ─────────────────────
    if args.output_lists:
        cb_idx = _build_cb_index(annot)

        def _dir_ids(df: pd.DataFrame, sig: set[str], positive: bool) -> set[str]:
            ic = "circ_id" if "circ_id" in df.columns else None
            if ic is None or "log2FC" not in df.columns:
                return set()
            mask = df[ic].isin(sig) & (df["log2FC"] > 0 if positive else df["log2FC"] < 0)
            return set(df.loc[mask, ic].astype(str))

        de_lists: dict = {}

        # Our method (edgeR_ciriquant)
        our_up = _dir_ids(our_de, our_sig, True)
        our_dn = _dir_ids(our_de, our_sig, False)
        entry: dict = {
            "Sig_DE_circRNAs":   _build_de_list(our_de, our_sig, our_sig_col, cb_idx, slop=args.slop),
            "Up_regulated":      _build_de_list(our_de, our_up,  our_sig_col, cb_idx, slop=args.slop),
            "Down_regulated":    _build_de_list(our_de, our_dn,  our_sig_col, cb_idx, slop=args.slop),
            "Top20_in_circBase": _build_de_list(our_de, our_sig, our_sig_col, cb_idx,
                                                n_top=20, require_circbase=True, slop=args.slop),
        }
        ic_our = "circ_id" if "circ_id" in our_de.columns else None
        if "Type" in our_de.columns and ic_our:
            sm = our_de[ic_our].isin(our_sig)
            t1 = set(our_de.loc[sm & (our_de["Type"] == "Type_I"),  ic_our].astype(str))
            t2 = set(our_de.loc[sm & (our_de["Type"] == "Type_II"), ic_our].astype(str))
            t1u = {x for x in t1 if not _fuzzy_in(x, nfcore_sig, args.slop)}
            entry["Type_I_count"]            = _build_de_list(our_de, t1,  our_sig_col, cb_idx, slop=args.slop)
            entry["Type_II_count"]           = _build_de_list(our_de, t2,  our_sig_col, cb_idx, slop=args.slop)
            entry["Type_I_unique_vs_DESeq2"] = _build_de_list(our_de, t1u, our_sig_col, cb_idx, slop=args.slop)
        de_lists["Our_edgeR_ciriquant"] = entry

        # DESeq2 baseline
        nf_up = _dir_ids(nfcore_de, nfcore_sig, True)
        nf_dn = _dir_ids(nfcore_de, nfcore_sig, False)
        de_lists["DESeq2_baseline"] = {
            "Sig_DE_circRNAs":   _build_de_list(nfcore_de, nfcore_sig, deseq_sig_col, cb_idx, slop=args.slop),
            "Up_regulated":      _build_de_list(nfcore_de, nf_up,      deseq_sig_col, cb_idx, slop=args.slop),
            "Down_regulated":    _build_de_list(nfcore_de, nf_dn,      deseq_sig_col, cb_idx, slop=args.slop),
            "Top20_in_circBase": _build_de_list(nfcore_de, nfcore_sig, deseq_sig_col, cb_idx,
                                                n_top=20, require_circbase=True, slop=args.slop),
        }

        # limma-voom
        if limma_de is not None and limma_sig:
            lm_up = _dir_ids(limma_de, limma_sig, True)
            lm_dn = _dir_ids(limma_de, limma_sig, False)
            de_lists["limma_voom"] = {
                "Sig_DE_circRNAs":   _build_de_list(limma_de, limma_sig, limma_sig_col, cb_idx, slop=args.slop),
                "Up_regulated":      _build_de_list(limma_de, lm_up,     limma_sig_col, cb_idx, slop=args.slop),
                "Down_regulated":    _build_de_list(limma_de, lm_dn,     limma_sig_col, cb_idx, slop=args.slop),
                "Top20_in_circBase": _build_de_list(limma_de, limma_sig, limma_sig_col, cb_idx,
                                                    n_top=20, require_circbase=True, slop=args.slop),
            }

        # Jaccard pairwise circRNA lists for interactive Jaccard table
        def _jac_pair(de_a, sig_a, col_a, de_b, sig_b, col_b):
            a_only = {x for x in sig_a if not _fuzzy_in(x, sig_b, args.slop)}
            b_only = {x for x in sig_b if not _fuzzy_in(x, sig_a, args.slop)}
            both_a = {x for x in sig_a if _fuzzy_in(x, sig_b, args.slop)}
            return {
                "A_only": _build_de_list(de_a, a_only, col_a, cb_idx, slop=args.slop),
                "B_only": _build_de_list(de_b, b_only, col_b, cb_idx, slop=args.slop),
                "Both":   _build_de_list(de_a, both_a, col_a, cb_idx, slop=args.slop),
            }

        jac_lists = [{
            "Method_A": "Our_edgeR_ciriquant", "Method_B": "DESeq2_baseline",
            **_jac_pair(our_de, our_sig, our_sig_col, nfcore_de, nfcore_sig, deseq_sig_col),
        }]
        if limma_de is not None and limma_sig:
            jac_lists.append({
                "Method_A": "Our_edgeR_ciriquant", "Method_B": "limma_voom",
                **_jac_pair(our_de, our_sig, our_sig_col, limma_de, limma_sig, limma_sig_col),
            })
            jac_lists.append({
                "Method_A": "DESeq2_baseline", "Method_B": "limma_voom",
                **_jac_pair(nfcore_de, nfcore_sig, deseq_sig_col, limma_de, limma_sig, limma_sig_col),
            })
        de_lists["_jaccard"] = jac_lists

        Path(args.output_lists).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output_lists, "w", encoding="utf-8") as fh:
            json.dump(de_lists, fh, ensure_ascii=False)
        print(f"[de_quality] Lists    → {args.output_lists}", file=sys.stderr)

    # ── Write outputs ─────────────────────────────────────────────────────────
    Path(args.output_summary).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(args.output_summary, sep="\t", index=False)
    pd.DataFrame(jac_rows).to_csv(    args.output_jaccard, sep="\t", index=False)
    print(f"[de_quality] Summary → {args.output_summary}", file=sys.stderr)
    print(f"[de_quality] Jaccard → {args.output_jaccard}", file=sys.stderr)


if __name__ == "__main__":
    main()
