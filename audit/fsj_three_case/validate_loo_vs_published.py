"""Gate: the LOO full-data (NONE) runs must reproduce the published DE tables before any LOO
output is used. Compares, per dataset and per method, log2FC and pvalue over all shared circ_ids.

usage (on the server, where the published results live):
  python3 validate_loo_vs_published.py <loo_root> <GSE> [<GSE> ...]
where <loo_root>/<GSE>/NONE/{edger,deseq2,limma}.tsv exists and the published tables are
~/<GSE>_results/de/{de_results.tsv, de_results_deseq2.tsv, de_results_limma.tsv}.

exit status 1 if any dataset/method differs by more than TOL in log2FC or pvalue, or if the
set of circ_ids differs. Nothing is written silently: every line is printed.
"""
import os
import sys
import pandas as pd

TOL = 1e-6
METHODS = [("edger.tsv", "de_results.tsv", "edgeR"),
           ("deseq2.tsv", "de_results_deseq2.tsv", "DESeq2"),
           ("limma.tsv", "de_results_limma.tsv", "limma")]

root = sys.argv[1]
gses = sys.argv[2:]
failed = []
for g in gses:
    for loo_f, pub_f, label in METHODS:
        loo_p = os.path.join(root, g, "NONE", loo_f)
        pub_p = os.path.join(os.path.expanduser("~"), g + "_results", "de", pub_f)
        if not os.path.exists(loo_p) or not os.path.exists(pub_p):
            print(f"{g:12s} {label:8s} MISSING ({loo_p} or {pub_p})"); failed.append((g, label)); continue
        loo = pd.read_csv(loo_p, sep="\t")
        pub = pd.read_csv(pub_p, sep="\t")
        if "log2FoldChange" in pub.columns:
            pub = pub.rename(columns={"log2FoldChange": "log2FC"})
        loo["circ_id"] = loo["circ_id"].astype(str); pub["circ_id"] = pub["circ_id"].astype(str)
        same_ids = set(loo.circ_id) == set(pub.circ_id)
        m = loo.merge(pub[["circ_id", "log2FC", "pvalue"]], on="circ_id", suffixes=("_loo", "_pub"))
        d_lfc = (m.log2FC_loo - m.log2FC_pub).abs().max() if len(m) else float("nan")
        d_p = (m.pvalue_loo - m.pvalue_pub).abs().max() if len(m) else float("nan")
        ok = same_ids and len(m) and d_lfc <= TOL and d_p <= TOL
        print(f"{g:12s} {label:8s} rows_loo={len(loo):6d} rows_pub={len(pub):6d} same_ids={same_ids} "
              f"max|dlog2FC|={d_lfc:.3g} max|dp|={d_p:.3g} {'OK' if ok else 'MISMATCH'}")
        if not ok:
            failed.append((g, label))
print("RESULT:", "all OK" if not failed else f"{len(failed)} MISMATCH/MISSING: {failed}")
sys.exit(0 if not failed else 1)
