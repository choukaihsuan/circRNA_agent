"""Parity check: biomarker scores shown in the PUBLISHED HTML reports vs the candidate files.

The published report.html embeds, per DE method, a bm_table computed inside generate_report.py
(a second scoring path). This script extracts that table and compares, for every circRNA in it,
the displayed biomarker_score / n_mirna / n_rbp with the rank_biomarkers candidate file for the
same method. Any difference means the two paths had already diverged for that dataset.

usage (server): python3 report_vs_candidate_parity.py <GSE> [<GSE> ...]
reads ~/<GSE>_results/report.html and ~/<GSE>_results/de/biomarker_candidates_v2_<method>.tsv
(v2 = current scoring; the v1 file is compared too when present).
"""
import json
import os
import re
import sys

import pandas as pd

METHOD_KEYS = {"edgeR_ciriquant": "edgeR", "deseq2": "DESeq2", "limma": "limma"}


def extract_all_de_methods(html: str) -> dict:
    m = re.search(r"ALL_DE_METHODS\s*=\s*", html)
    if not m:
        raise SystemExit("ALL_DE_METHODS not found in report")
    obj, _ = json.JSONDecoder().raw_decode(html, m.end())
    return obj


for gse in sys.argv[1:]:
    base = os.path.expanduser(f"~/{gse}_results")
    html = open(os.path.join(base, "report.html"), encoding="utf-8").read()
    allm = extract_all_de_methods(html)
    print(f"== {gse}: report methods {sorted(allm)}")
    for rkey, meth in METHOD_KEYS.items():
        if rkey not in allm or "bm_table" not in allm[rkey]:
            print(f"   {meth}: no bm_table in report"); continue
        bm = allm[rkey]["bm_table"]
        cols = bm["cols"]
        rows = pd.DataFrame(bm["rows"], columns=cols)
        rows["circ_id"] = rows["circ_id"].astype(str)
        for label, cand_p in [("v2", os.path.join(base, "de", f"biomarker_candidates_v2_{meth}.tsv")),
                              ("v1", os.path.join(base, "de", "biomarker_candidates.tsv") if meth == "edgeR" else None)]:
            if not cand_p or not os.path.exists(cand_p):
                continue
            cand = pd.read_csv(cand_p, sep="\t")
            cand["circ_id"] = cand["circ_id"].astype(str)
            m2 = rows.merge(cand[["circ_id", "biomarker_score", "n_mirna", "n_rbp"]],
                            on="circ_id", how="left", suffixes=("_report", "_file"))
            missing = int(m2["biomarker_score_file"].isna().sum())
            d_score = (m2.biomarker_score_report.astype(float) - m2.biomarker_score_file.astype(float)).abs()
            d_m = (m2.n_mirna_report.astype(float) - m2.n_mirna_file.astype(float)).abs()
            d_r = (m2.n_rbp_report.astype(float) - m2.n_rbp_file.astype(float)).abs()
            n_diff = int((d_score > 1e-3).sum())
            print(f"   {meth:7s} vs {label}: report_rows={len(rows):5d} not_in_file={missing:3d} "
                  f"score_diff_rows={n_diff:5d} max|dscore|={d_score.max():.4g} "
                  f"max|dn_mirna|={d_m.max():.4g} max|dn_rbp|={d_r.max():.4g}")
