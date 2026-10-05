"""Terminal regeneration for the 15 analysed datasets: rank_biomarkers_v2 (3 methods) and generate_report.

No Snakemake DAG. Every argument mirrors workflow/rules/de.smk for the same project config; the
LOO directory, cross-method minimum and stability ratio come from the datasets' own status files
(the values the published candidates were produced with). --no-run-params is passed so the
run_params.yaml file is not created here (that is item B5 and is not part of this step).

usage: regen.py <repo_dir> <GSE> [<GSE> ...]
"""
import json
import os
import shlex
import subprocess
import sys
import types

import yaml

PY = os.environ.get("PYTHON", "python")
repo = os.path.abspath(sys.argv[1])
os.chdir(repo)

PRIMARY = {"edgeR": "de_results.tsv", "DESeq2": "de_results_deseq2.tsv", "limma": "de_results_limma.tsv"}
LOG = os.path.expanduser("~/item8_regen/logs")
os.makedirs(LOG, exist_ok=True)


def run(cmd, log_name):
    with open(os.path.join(LOG, log_name), "w") as lf:
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0:
        raise SystemExit(f"FAILED rc={rc}: {' '.join(shlex.quote(c) for c in cmd)}")


for gse in sys.argv[2:]:
    cfg = yaml.safe_load(open(f"config/projects/{gse}.yaml"))
    res = cfg["results_dir"]
    de = cfg["de"]
    st = json.load(open(f"{res}/de/biomarker_candidates_v2_edgeR.status.json"))
    loo_dir, cmm, smr = st["loo_dir"], st["cross_method_min"], st["stability_min_ratio"]
    for m in ["edgeR", "DESeq2", "limma"]:
        s = json.load(open(f"{res}/de/biomarker_candidates_v2_{m}.status.json"))
        assert s["loo_dir"] == loo_dir and s["cross_method_min"] == cmm and s["stability_min_ratio"] == smr, \
            f"{gse} {m}: status files disagree"
    # the rule's inputs (de.smk, rank_biomarkers_v2)
    cmd = [PY, "scripts/rank_biomarkers.py",
           "--de", f"{res}/de/{PRIMARY['edgeR']}",
           "--primary-method", "edgeR"]  # replaced per method below
    for m in ["edgeR", "DESeq2", "limma"]:
        cmd = [PY, "scripts/rank_biomarkers.py",
               "--de", f"{res}/de/{PRIMARY[m]}", "--primary-method", m,
               "--de-edger", f"{res}/de/de_results.tsv",
               "--de-deseq2", f"{res}/de/de_results_deseq2.tsv",
               "--de-limma", f"{res}/de/de_results_limma.tsv",
               "--annot", f"{res}/circRNA/circbase_annotated.tsv",
               "--isoforms", f"{res}/circRNA/isoform_groups.tsv",
               "--interactions", f"{res}/de/interactions.json",
               "--output", f"{res}/de/biomarker_candidates_v2_{m}.tsv",
               "--status-json", f"{res}/de/biomarker_candidates_v2_{m}.status.json",
               "--fdr", str(de["fdr_cutoff"]), "--lfc", str(de["log2fc_cutoff"]),
               "--de-sig-by", de.get("de_sig_by", "pvalue"),
               "--loo-dir", loo_dir, "--cross-method-min", str(cmm),
               "--stability-min-ratio", str(smr), "--no-run-params"]
        run(cmd, f"{gse}_rank_{m}.log")
        print("rank ok", gse, m, flush=True)

    # generate_report: mock Snakemake object, same inputs and params as de.smk generate_report
    gr_in = {
        "de": f"{res}/de/de_results.tsv",
        "de_edger": f"{res}/de/de_results_edgeR_ciriquant.tsv",
        "de_deseq": f"{res}/de/de_results_deseq2.tsv",
        "de_limma": f"{res}/de/de_results_limma.tsv",
        "biomarkers_edger": f"{res}/de/biomarker_candidates_v2_edgeR.tsv",
        "biomarkers_deseq2": f"{res}/de/biomarker_candidates_v2_DESeq2.tsv",
        "biomarkers_limma": f"{res}/de/biomarker_candidates_v2_limma.tsv",
        "biomarker_status": [f"{res}/de/biomarker_candidates_v2_{m}.status.json" for m in ["edgeR", "DESeq2", "limma"]],
        "matrix": f"{res}/circRNA/count_matrix.tsv",
        "volcano": f"{res}/plots/volcano.pdf",
        "heatmap": f"{res}/plots/heatmap.pdf",
        "pca": f"{res}/plots/pca.pdf",
        "multiqc": f"{res}/qc/multiqc_report.html",
        "switching": f"{res}/de/isoform_switching.tsv",
        "groups": cfg["groups"],
        "isoform_groups": f"{res}/circRNA/isoform_groups.tsv",
        "circbase_annot": f"{res}/circRNA/circbase_annotated.tsv",
        "interactions": f"{res}/de/interactions.json",
    }
    de_method = cfg.get("de", {}).get("method", "edgeR_ciriquant")
    gr_params = {
        "project_id": cfg["project_id"], "fdr": de["fdr_cutoff"], "lfc": de["log2fc_cutoff"],
        "de_method": de_method, "de_sig_by": de.get("de_sig_by", "pvalue"),
        "tumor_label": de["tumor_label"], "normal_label": de["normal_label"],
        "heatmap_top_n": de.get("heatmap_top_n", 10), "study_title": cfg.get("study_title", ""),
    }
    sn = types.SimpleNamespace()
    sn.input = types.SimpleNamespace(**gr_in)
    sn.params = types.SimpleNamespace(**gr_params)
    sn.output = [f"{res}/report.html"]
    sn.log = [os.path.join(LOG, f"{gse}_generate_report.log")]
    code = open("scripts/generate_report.py", encoding="utf-8").read()
    with open(sn.log[0], "w") as lf:
        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = lf
        try:
            exec(compile(code, "scripts/generate_report.py", "exec"), {"snakemake": sn, "__name__": "__main__"})
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr
    print("report ok", gse, flush=True)
print("REGEN_DONE", flush=True)
