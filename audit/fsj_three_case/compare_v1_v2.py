"""Old (v1 biomarker_candidates.tsv) vs new (v2 biomarker_candidates_v2.tsv) top-10 comparison.

For each dataset: v1 top-10 by score (the published ranking), v2 top-10 under the default
view (rank_default <= 10: three-method AND effect-stable), overlap, and for every v1 top-10
member the reason it would not appear in the v2 default view.

usage: python3 compare_v1_v2.py <dir with <GSE>/biomarker_candidates.tsv and <GSE>/biomarker_candidates_v2.tsv> <out.tsv>
"""
import os
import sys
import pandas as pd

root, out = sys.argv[1], sys.argv[2]
GSES = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
        "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289", "SRP156355"]


def reason(row):
    if not bool(row["three_method_sig"]):
        return "not_three_method(n_sig=%s)" % row["n_sig_methods"]
    bad = [m for m in ("edgeR", "DESeq2", "limma") if row.get("stable_%s" % m) is False]
    if bad:
        return "unstable(" + ",".join(bad) + ")"
    return "pass"


rows, detail = [], []
for g in GSES:
    v1p = os.path.join(root, g, "biomarker_candidates.tsv")
    v2p = os.path.join(root, g, "biomarker_candidates_v2.tsv")
    if not (os.path.exists(v1p) and os.path.exists(v2p)):
        rows.append(dict(dataset=g, error="missing v1 or v2 file")); continue
    v1 = pd.read_csv(v1p, sep="\t")
    v2 = pd.read_csv(v2p, sep="\t")
    v1 = v1.sort_values("rank").head(10)
    v1_ids = v1["circ_id"].astype(str).tolist()
    v2["circ_id"] = v2["circ_id"].astype(str)
    v2_top = v2[v2["rank_default"].notna() & (pd.to_numeric(v2["rank_default"], errors="coerce") <= 10)]
    v2_top = v2_top.sort_values("rank_default")
    v2_ids = v2_top["circ_id"].tolist()
    overlap = [c for c in v1_ids if c in v2_ids]
    reasons = {}
    for c in v1_ids:
        r = v2.loc[v2.circ_id == c]
        if r.empty:
            reasons.setdefault("absent_in_v2", 0); reasons["absent_in_v2"] += 1; continue
        k = reason(r.iloc[0])
        reasons[k] = reasons.get(k, 0) + 1
        detail.append(dict(dataset=g, circ_id=c, circbase_gene=r.iloc[0].get("circbase_gene", ""),
                           reason=k, v1_rank=int(v1.loc[v1.circ_id.astype(str) == c, "rank"].iloc[0])))
    n_pass = int(v2["pass_default"].astype(str).eq("True").sum())
    rows.append(dict(dataset=g, v1_candidates=len(pd.read_csv(v1p, sep="\t")), v2_candidates=len(v2),
                     v2_three_method=int(v2["three_method_sig"].astype(str).eq("True").sum()),
                     v2_pass_default=n_pass, v1_top10_in_v2_top10=len(overlap),
                     v1_top10_reasons="; ".join("%s=%d" % kv for kv in sorted(reasons.items())),
                     v2_top10=";".join(v2_ids), v1_top10=";".join(v1_ids)))

summary = pd.DataFrame(rows)
summary.to_csv(out, sep="\t", index=False)
pd.DataFrame(detail).to_csv(out.replace(".tsv", "_top10_detail.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print(summary.drop(columns=[c for c in ("v1_top10", "v2_top10") if c in summary.columns]).to_string(index=False))
