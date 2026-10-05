"""Column-level check of regenerated v2 candidate files against the archived pre-regeneration copies.

Expected: same row set, same columns, same values in every column except `rank`; `rank` may differ
only inside groups of exactly equal biomarker_score. rank_default and pass_default must be identical.
"""
import os
import subprocess
import sys
import tarfile

import pandas as pd

ARCH = os.path.expanduser("~/archive_published_20261005")
HOME = os.path.expanduser("~")
GSES = ["GSE108735", "GSE113230", "GSE121842", "GSE130078", "GSE133998", "GSE136569", "GSE143797",
        "GSE148036", "GSE171011", "GSE221107", "GSE248612", "GSE77509", "GSE97239", "PRJNA553289", "SRP156355"]
totals = dict(files=0, rank_changed=0, other_col_changed=0, rank_default_changed=0, unexplained=0)
for g in GSES:
    tf = tarfile.open(os.path.join(ARCH, f"{g}_results.tar.gz"))
    for m in ["edgeR", "DESeq2", "limma"]:
        member = f"{g}_results/de/biomarker_candidates_v2_{m}.tsv"
        old = pd.read_csv(tf.extractfile(member), sep="\t", dtype=str, keep_default_na=False)
        new = pd.read_csv(os.path.join(HOME, f"{g}_results/de/biomarker_candidates_v2_{m}.tsv"),
                          sep="\t", dtype=str, keep_default_na=False)
        totals["files"] += 1
        assert list(old.columns) == list(new.columns), f"{g} {m}: columns differ"
        assert len(old) == len(new), f"{g} {m}: row count differs"
        o = old.set_index("circ_id")
        n = new.set_index("circ_id").loc[o.index]
        assert set(old.circ_id) == set(new.circ_id), f"{g} {m}: circ_id set differs"
        other = [c for c in o.columns if c != "rank" and not (o[c] == n[c]).all()]
        if other:
            totals["other_col_changed"] += 1
            print(f"UNEXPLAINED {g} {m}: columns changed {other}")
        rk = o["rank"] != n["rank"]
        totals["rank_changed"] += int(rk.sum())
        for col in ["rank_default", "pass_default"]:
            if not (o[col] == n[col]).all():
                totals["rank_default_changed"] += 1
                print(f"UNEXPLAINED {g} {m}: {col} changed")
        # every rank change must sit in a group of exactly equal biomarker_score
        s = o["biomarker_score"].astype(float)
        tied = s.map(s.value_counts()) > 1
        bad = int((rk & ~tied).sum())
        totals["unexplained"] += bad
        if bad:
            print(f"UNEXPLAINED {g} {m}: {bad} rank changes outside exact ties")
        print(f"{g:11s} {m:7s} rows={len(o):5d} rank_changed={int(rk.sum()):4d} "
              f"rank_default_identical={bool((o.rank_default == n.rank_default).all())}")
print("TOTAL", totals)
