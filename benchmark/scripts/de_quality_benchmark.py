"""
de_quality_benchmark.py – DE 分析品質比較

比較三種方法對 GSE113230（肝癌 6 samples）的 DE 分析結果：
  1. Our method      : edgeR_ciriquant (BSJ/FSJ ratio + Type I/II 分類)
  2. nf-core sim     : DESeq2 on consensus circRNA BSJ counts (no FSJ offset)
  3. sponging sim    : DESeq2 on DCC-only circRNA BSJ counts (optional)

計算指標：
  - 各方法顯著 DE circRNA 數量
  - Jaccard similarity（兩兩重疊）
  - Type I circRNA 中，只有我們方法找到的比例（nf-core 未見）
  - Top 20 DE circRNA 中的 circBase 已知比例

Outputs:
  --output-summary   results/benchmark/de_quality_summary.tsv
  --output-jaccard   results/benchmark/de_jaccard.tsv
"""

from __future__ import annotations

import argparse
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


def _sig_ids(df: pd.DataFrame, fdr: float, lfc: float) -> set[str]:
    """Return significant circRNA IDs."""
    if "padj" not in df.columns or "log2FC" not in df.columns:
        return set()
    mask = (df["padj"] < fdr) & (df["log2FC"].abs() > lfc)
    if "circ_id" in df.columns:
        return set(df.loc[mask, "circ_id"].dropna().astype(str))
    return set(df.loc[mask].index.astype(str))


def _circbase_hits(
    de: pd.DataFrame,
    annot: pd.DataFrame | None,
    n_top: int = 20,
    slop: int = 10,
) -> int | None:
    """Count known circBase entries among top n_top DE circRNAs."""
    if annot is None or annot.empty:
        return None
    sorted_de = de.sort_values("padj").head(n_top)
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

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare DE analysis quality across pipeline simulations"
    )
    parser.add_argument("--our-de",       required=True,
                        help="edgeR_ciriquant DE results TSV")
    parser.add_argument("--nfcore-de",    required=True,
                        help="DESeq2 DE results (nf-core sim)")
    parser.add_argument("--sponging-de",  default=None,
                        help="DESeq2 DE results (sponging sim; optional)")
    parser.add_argument("--circbase-annot", default=None,
                        help="circbase_annotated.tsv from annotate_circbase.py")
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-jaccard", required=True)
    parser.add_argument("--fdr",  type=float, default=0.05)
    parser.add_argument("--lfc",  type=float, default=1.0)
    parser.add_argument("--slop", type=int,   default=10)
    args = parser.parse_args()

    # ── Load DE results ───────────────────────────────────────────────────────
    our_de     = _normalise_cols(pd.read_csv(args.our_de, sep="\t"))
    nfcore_de  = _normalise_cols(pd.read_csv(args.nfcore_de, sep="\t"))
    sponging_de: pd.DataFrame | None = None
    if args.sponging_de and Path(args.sponging_de).exists():
        sponging_de = _normalise_cols(pd.read_csv(args.sponging_de, sep="\t"))

    annot: pd.DataFrame | None = None
    if args.circbase_annot and Path(args.circbase_annot).exists():
        annot = pd.read_csv(args.circbase_annot, sep="\t")

    # ── Significant circRNA sets ──────────────────────────────────────────────
    our_sig     = _sig_ids(our_de,    args.fdr, args.lfc)
    nfcore_sig  = _sig_ids(nfcore_de, args.fdr, args.lfc)
    sponge_sig  = _sig_ids(sponging_de, args.fdr, args.lfc) if sponging_de is not None else set()

    print(
        f"[de_quality] Significant (FDR<{args.fdr}, |lFC|>{args.lfc}): "
        f"ours={len(our_sig)}  nfcore={len(nfcore_sig)}  sponging={len(sponge_sig)}",
        file=sys.stderr,
    )

    # ── Type I / II counts ────────────────────────────────────────────────────
    n_type1 = n_type2 = n_type1_unique = None
    if "Type" in our_de.columns:
        id_col = "circ_id" if "circ_id" in our_de.columns else None
        sig_de  = our_de[
            our_de[id_col].isin(our_sig) if id_col
            else our_de.index.isin(our_sig)
        ]
        n_type1 = int((sig_de["Type"] == "Type_I").sum())
        n_type2 = int((sig_de["Type"] == "Type_II").sum())

        # Type I unique: not found in nf-core significant set
        type1_ids = set(
            sig_de.loc[sig_de["Type"] == "Type_I", id_col].astype(str).tolist()
            if id_col else
            sig_de[sig_de["Type"] == "Type_I"].index.astype(str).tolist()
        )
        n_type1_unique = sum(
            1 for cid in type1_ids
            if not _fuzzy_in(cid, nfcore_sig, args.slop)
        )
        print(
            f"[de_quality] Type I={n_type1}  Type II={n_type2}  "
            f"Type I unique vs nfcore={n_type1_unique}",
            file=sys.stderr,
        )

    # ── circBase hits ─────────────────────────────────────────────────────────
    cb_our     = _circbase_hits(our_de,     annot, 20, args.slop)
    cb_nfcore  = _circbase_hits(nfcore_de,  annot, 20, args.slop)
    cb_sponge  = _circbase_hits(sponging_de, annot, 20, args.slop) \
                 if sponging_de is not None else None

    # ── Summary table ─────────────────────────────────────────────────────────
    summary_rows = [
        {
            "Method":                  "Our_edgeR_ciriquant",
            "Total_input_circRNAs":    len(our_de),
            "Sig_DE_circRNAs":         len(our_sig),
            "Up_regulated":            int(
                (our_de.loc[
                    our_de["circ_id"].isin(our_sig)
                    if "circ_id" in our_de.columns
                    else our_de.index.isin(our_sig),
                    "log2FC"
                ] > 0).sum()
            ) if "log2FC" in our_de.columns else None,
            "Down_regulated":          int(
                (our_de.loc[
                    our_de["circ_id"].isin(our_sig)
                    if "circ_id" in our_de.columns
                    else our_de.index.isin(our_sig),
                    "log2FC"
                ] < 0).sum()
            ) if "log2FC" in our_de.columns else None,
            "Type_I_count":            n_type1,
            "Type_II_count":           n_type2,
            "Type_I_unique_vs_nfcore": n_type1_unique,
            "Top20_in_circBase":       cb_our,
            "DE_method":               "edgeR_ciriquant (BSJ/FSJ ratio + offset)",
        },
        {
            "Method":                  "nfcore_DESeq2",
            "Total_input_circRNAs":    len(nfcore_de),
            "Sig_DE_circRNAs":         len(nfcore_sig),
            "Up_regulated":            int(
                (nfcore_de.loc[
                    nfcore_de["circ_id"].isin(nfcore_sig)
                    if "circ_id" in nfcore_de.columns
                    else nfcore_de.index.isin(nfcore_sig),
                    "log2FC"
                ] > 0).sum()
            ) if "log2FC" in nfcore_de.columns else None,
            "Down_regulated":          int(
                (nfcore_de.loc[
                    nfcore_de["circ_id"].isin(nfcore_sig)
                    if "circ_id" in nfcore_de.columns
                    else nfcore_de.index.isin(nfcore_sig),
                    "log2FC"
                ] < 0).sum()
            ) if "log2FC" in nfcore_de.columns else None,
            "Type_I_count":            "N/A",
            "Type_II_count":           "N/A",
            "Type_I_unique_vs_nfcore": "N/A",
            "Top20_in_circBase":       cb_nfcore,
            "DE_method":               "DESeq2 (BSJ counts only, simulated)",
        },
    ]

    if sponging_de is not None:
        summary_rows.append({
            "Method":                  "sponging_DESeq2",
            "Total_input_circRNAs":    len(sponging_de),
            "Sig_DE_circRNAs":         len(sponge_sig),
            "Up_regulated":            int(
                (sponging_de.loc[
                    sponging_de["circ_id"].isin(sponge_sig)
                    if "circ_id" in sponging_de.columns
                    else sponging_de.index.isin(sponge_sig),
                    "log2FC"
                ] > 0).sum()
            ) if "log2FC" in sponging_de.columns else None,
            "Down_regulated":          int(
                (sponging_de.loc[
                    sponging_de["circ_id"].isin(sponge_sig)
                    if "circ_id" in sponging_de.columns
                    else sponging_de.index.isin(sponge_sig),
                    "log2FC"
                ] < 0).sum()
            ) if "log2FC" in sponging_de.columns else None,
            "Type_I_count":            "N/A",
            "Type_II_count":           "N/A",
            "Type_I_unique_vs_nfcore": "N/A",
            "Top20_in_circBase":       cb_sponge,
            "DE_method":               "DESeq2 (DCC-only BSJ, simulated)",
        })

    # ── Jaccard pairwise table ────────────────────────────────────────────────
    jac_rows = [
        {
            "Method_A":  "Our_edgeR_ciriquant",
            "Method_B":  "nfcore_DESeq2",
            "Jaccard":   _jaccard(our_sig, nfcore_sig, args.slop),
            "A_only":    sum(1 for x in our_sig if not _fuzzy_in(x, nfcore_sig, args.slop)),
            "B_only":    sum(1 for x in nfcore_sig if not _fuzzy_in(x, our_sig, args.slop)),
            "Both":      sum(1 for x in our_sig if _fuzzy_in(x, nfcore_sig, args.slop)),
        },
    ]
    if sponge_sig:
        jac_rows.append({
            "Method_A":  "Our_edgeR_ciriquant",
            "Method_B":  "sponging_DESeq2",
            "Jaccard":   _jaccard(our_sig, sponge_sig, args.slop),
            "A_only":    sum(1 for x in our_sig    if not _fuzzy_in(x, sponge_sig, args.slop)),
            "B_only":    sum(1 for x in sponge_sig if not _fuzzy_in(x, our_sig,    args.slop)),
            "Both":      sum(1 for x in our_sig    if     _fuzzy_in(x, sponge_sig, args.slop)),
        })
        jac_rows.append({
            "Method_A":  "nfcore_DESeq2",
            "Method_B":  "sponging_DESeq2",
            "Jaccard":   _jaccard(nfcore_sig, sponge_sig, args.slop),
            "A_only":    sum(1 for x in nfcore_sig if not _fuzzy_in(x, sponge_sig, args.slop)),
            "B_only":    sum(1 for x in sponge_sig if not _fuzzy_in(x, nfcore_sig, args.slop)),
            "Both":      sum(1 for x in nfcore_sig if     _fuzzy_in(x, sponge_sig, args.slop)),
        })

    # ── Write outputs ─────────────────────────────────────────────────────────
    Path(args.output_summary).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows).to_csv(args.output_summary, sep="\t", index=False)
    pd.DataFrame(jac_rows).to_csv(    args.output_jaccard, sep="\t", index=False)
    print(f"[de_quality] Summary → {args.output_summary}", file=sys.stderr)
    print(f"[de_quality] Jaccard → {args.output_jaccard}", file=sys.stderr)


if __name__ == "__main__":
    main()
