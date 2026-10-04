"""
partner_comparison.py – CIRIquant 固定為主 caller，比較候選搭檔的佐證資訊量。

對 GSE55872 SRR444655 的現成偵測輸出（不重跑任何比對）執行 consensus_filter.py，
再對 RNase R ground truth 評估。獨立於 accuracy_summary.tsv，不修改任何既有 benchmark 輸出。

  prec_corroborated    = 被搭檔佐證的 CIRIquant 判定之 precision（= 共識輸出）
  prec_uncorroborated  = 未被佐證的 CIRIquant 判定之 precision（CIRIquant 單獨 − 共識）
  lift                 = prec_corroborated − CIRIquant 單獨 precision
  discrimination       = prec_corroborated − prec_uncorroborated

評估為 ground-truth 導向（與 accuracy_benchmark.evaluate 相同）：以 GT 列是否落在偵測集合
slop 範圍內計 TP/FP/FN/TN；uncorroborated 集合 = GT 列被 CIRIquant 偵測、但未被共識集合偵測者。
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from accuracy_benchmark import _coord_match_in_set, _load_bed, _prf  # noqa: E402


def run_consensus(cf: str, out_prefix: Path, inputs: dict, min_tools: int, slop: int, min_bsj: int):
    bed, summ = f"{out_prefix}.bed", f"{out_prefix}.summary.tsv"
    cmd = [sys.executable, cf, "--output", bed, "--summary", summ,
           "--min-tools", str(min_tools), "--slop", str(slop),
           "--min-bsj", str(min_bsj), "--max-junction-ratio", "999"]
    for flag, path in inputs.items():
        cmd += [flag, path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    out_prefix.with_suffix(".log").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"consensus_filter failed ({out_prefix.name}):\n{r.stderr[-800:]}")
    return _load_bed(bed)


def score(ids_a, gt, slop, subtract=None):
    """TP/FP/FN/TN over GT rows. If subtract is given, rows detected by ids_a but NOT by subtract."""
    tp = fp = fn = tn = 0
    for _, row in gt.iterrows():
        det = _coord_match_in_set(row["circ_id"], ids_a, slop)
        if subtract is not None:
            det = det and not _coord_match_in_set(row["circ_id"], subtract, slop)
        lbl = int(row["is_true"])
        if lbl == 1 and det: tp += 1
        elif lbl == 0 and det: fp += 1
        elif lbl == 1: fn += 1
        else: tn += 1
    return tp, fp, fn, tn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--consensus-filter", required=True)
    ap.add_argument("--sample-dir", required=True, help="…/GSE55872_results/circRNA/SRR444655")
    ap.add_argument("--ground-truth", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--slop", type=int, default=10)
    ap.add_argument("--min-bsj", type=int, default=2)
    a = ap.parse_args()

    sd = Path(a.sample_dir)
    wd = Path(a.workdir); wd.mkdir(parents=True, exist_ok=True)
    CIRI = {"--cirique": str(sd / "SRR444655.gtf")}
    DCC = {"--dcc": str(sd / "DCC/CircCoordinates")}
    CE2 = {"--circexplorer2": str(sd / "CIRCexplorer2/known_circ.txt")}
    FC = {"--find-circ": str(sd / "find_circ/splice_sites.bed")}  # parse_find_circ keeps CIRCULAR only

    gt = pd.read_csv(a.ground_truth, sep="\t")
    gt = gt[gt["is_true"].isin([0, 1])]

    def go(name, inputs, mt):
        return run_consensus(a.consensus_filter, wd / name, inputs, mt, a.slop, a.min_bsj)

    ciri = go("ciriquant_only", CIRI, 1)
    dcc_only = go("dcc_only", DCC, 1)
    partners = [("DCC", DCC, False), ("CIRCexplorer2", CE2, False), ("find_circ", FC, True)]

    def row(partner, ids, excl, base_prec=None, corrob=None):
        tp, fp, fn, tn = score(ids, gt, a.slop)
        p, r, f = _prf(tp, fp, fn)
        d = dict(partner=partner, detected=len(ids), TP=tp, FP=fp, FN=fn, TN=tn, precision=p,
                 recall=r, F1=f, specificity=round(tn / (tn + fp), 4) if tn + fp else 0.0,
                 prec_corroborated=None, prec_uncorroborated=None, lift=None,
                 discrimination=None, retained_calls=len(ids), excluded_by_criterion=excl)
        if corrob is not None:
            utp, ufp, _, _ = score(ciri, gt, a.slop, subtract=ids)
            up = round(utp / (utp + ufp), 4) if utp + ufp else None
            d.update(prec_corroborated=p, prec_uncorroborated=up,
                     lift=round(p - base_prec, 4),
                     discrimination=None if up is None else round(p - up, 4),
                     uncorrob_TP=utp, uncorrob_FP=ufp)
        return d

    rows = [row("CIRIquant_only", ciri, False)]
    base = rows[0]["precision"]
    for name, inp, excl in partners:
        ids = go(f"ciriquant_plus_{name}", {**CIRI, **inp}, 2)
        rows.append(row(f"CIRIquant+{name}", ids, excl, base_prec=base, corrob=True))
    rows.append(row("DCC_only", dcc_only, False))
    df = pd.DataFrame(rows)
    df.to_csv(a.output, sep="\t", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
