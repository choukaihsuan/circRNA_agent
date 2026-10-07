"""Figure 12 data: the 6 raw/normalized biomarker-score components (S,F,C,K,M,R) for the edgeR
significant set of each dataset, built by importing rank_biomarkers.py's own functions (not a
reimplementation) so the per-dataset, per-locus normalization matches the published pipeline
exactly. Validates against the published biomarker_score for GSE133998 before using the values.

usage: fig12_components.py <HOME>
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HOME = sys.argv[1]
sys.path.insert(0, os.path.join(HOME, "circRNA_agent", "scripts"))
import rank_biomarkers as rb  # noqa: E402

GSES = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
        "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289",
        "SRP156355"]
LETTER = {"sig": "S", "fc": "F", "conf": "C", "known": "K", "mirna": "M", "rbp": "R"}
SIG_CAP, FC_CAP = 10.0, 5.0


def build_sig_and_norm(G):
    """Exact replica of rank_biomarkers.main()'s sig-building + the raw/norm block of
    _compute_weighted_score (copied verbatim from that function so the two cannot diverge
    silently; validated below against the published biomarker_score)."""
    base = os.path.join(HOME, f"{G}_results")
    de = pd.read_csv(os.path.join(base, "de", "de_results_edgeR_ciriquant.tsv"), sep="\t")
    annot = pd.read_csv(os.path.join(base, "circRNA", "circbase_annotated.tsv"), sep="\t")
    annot_cols = [c for c in ["circ_id", "confidence_score", "in_circbase", "circbase_id", "circbase_gene"]
                  if c in annot.columns]
    merged = de.merge(annot[annot_cols], on="circ_id", how="left")
    merged["in_circbase"] = merged.get("in_circbase", pd.Series(0, index=merged.index)).fillna(0).astype(int)
    merged["confidence_score"] = pd.to_numeric(merged.get("confidence_score", 0), errors="coerce").fillna(0)

    sig_col, sig_thr = rb._resolve_sig(merged, "pvalue", 0.05)
    sig_mask = merged[sig_col].notna() & (merged[sig_col] < sig_thr) & (merged["log2FC"].abs() > 1.0)
    sig = merged.loc[sig_mask].copy()

    icounts = rb._load_interaction_counts(os.path.join(base, "de", "interactions.json"))
    queried_ids = set(icounts.keys())
    sig["interaction_queried"] = sig["circ_id"].isin(queried_ids)
    sig["n_mirna"] = sig["circ_id"].map(lambda x: icounts.get(x, (0, 0))[0])
    sig["n_rbp"] = sig["circ_id"].map(lambda x: icounts.get(x, (0, 0))[1])

    # the real function (used for validation and as the score of record)
    score = rb._compute_weighted_score(sig, sig_col, SIG_CAP, FC_CAP, None)

    # raw/norm block, copied verbatim from _compute_weighted_score for the per-component values
    raw = pd.DataFrame(index=sig.index)
    raw["sig"] = sig[sig_col].apply(lambda p: min(-np.log10(max(p, 1e-10)), SIG_CAP))
    raw["fc"] = sig["log2FC"].abs().clip(upper=FC_CAP)
    raw["conf"] = sig["confidence_score"].clip(lower=0)
    raw["known"] = sig["in_circbase"].astype(float)
    raw["mirna"] = sig["n_mirna"].where(sig["interaction_queried"])
    raw["rbp"] = sig["n_rbp"].where(sig["interaction_queried"])
    norm = pd.DataFrame(index=sig.index, columns=raw.columns, dtype=float)
    for c in raw.columns:
        if c == "known":
            norm[c] = raw[c]
            continue
        valid = raw[c].notna()
        if valid.any():
            norm.loc[valid, c] = rb._minmax(raw.loc[valid, c])

    return sig, score, norm


if __name__ == "__main__":
    # validation on GSE133998
    sig, score, norm = build_sig_and_norm("GSE133998")
    pub = pd.read_csv(os.path.join(HOME, "GSE133998_results", "de", "biomarker_candidates.tsv"), sep="\t")
    pub = pub.set_index("circ_id")["biomarker_score"]
    mine = score.copy(); mine.index = sig["circ_id"].values
    common = pub.index.intersection(mine.index)
    diff = (pub.loc[common] - mine.loc[common]).abs()
    print(f"VALIDATION GSE133998: n_sig={len(sig)} n_common={len(common)} max|diff|={diff.max():.6f} "
          f"(0 = exact match to the published pipeline)")
    if diff.max() > 1e-3:
        print("STOP: reconstruction does not match the published pipeline; not using these values.")
        sys.exit(2)

    # per-dataset medians for all 15
    print()
    print(f"{'dataset':12s} {'n_sig':>6s} " + " ".join(f"{LETTER[c]:>7s}" for c in ["sig","fc","conf","known","mirna","rbp"]) + "  top_component")
    rows = []
    for G in GSES:
        s, sc, nm = build_sig_and_norm(G)
        med = nm.median()
        top = med.idxmax()
        rows.append((G, len(s), med, top))
        print(f"{G:12s} {len(s):6d} " + " ".join(f"{med[c]:7.4f}" for c in ["sig","fc","conf","known","mirna","rbp"]) + f"  {LETTER[top]} ({top})")

    tops = [LETTER[r[3]] for r in rows]
    from collections import Counter
    print()
    print("top-component tally across 15 datasets:", dict(Counter(tops)))

    json.dump({G: {"n_sig": n, "median": med.to_dict(), "top": top} for G, n, med, top in rows},
              open(os.path.join(HOME, "item16", "fig12_medians.json"), "w"), indent=1)

    # save GSE133998's norm + score for the figure script
    sig133998, score133998, norm133998 = build_sig_and_norm("GSE133998")
    out = norm133998.copy()
    out["circ_id"] = sig133998["circ_id"].values
    out["biomarker_score"] = score133998.values
    out.to_csv(os.path.join(HOME, "item16", "fig12_gse133998_components.tsv"), sep="\t", index=False)
    print("\nwrote", os.path.join(HOME, "item16", "fig12_gse133998_components.tsv"))
