"""
Differential expression and reporting rules (requires sample_groups.csv).
Skipped automatically when HAS_GROUPS is False (see Snakefile).
"""


rule de_analysis:
    """circRNA DE analysis; method = edgeR_ciriquant / deseq2 / limma (config de.method)."""
    input:
        matrix         = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        fsj_matrix     = RESULTS_DIR + "/circRNA/fsj_count_matrix.tsv",
        groups         = config["groups"],
        circbase_annot = RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
    output:
        de       = RESULTS_DIR + "/de/de_results.tsv",
        de_edger = RESULTS_DIR + "/de/de_results_edgeR_ciriquant.tsv",
        de_deseq = RESULTS_DIR + "/de/de_results_deseq2.tsv",
        de_limma = RESULTS_DIR + "/de/de_results_limma.tsv",
        volcano  = RESULTS_DIR + "/plots/volcano.pdf",
        heatmap  = RESULTS_DIR + "/plots/heatmap.pdf",
        pca      = RESULTS_DIR + "/plots/pca.pdf",
    params:
        de_method            = DE_METHOD,
        fdr                  = config["de"]["fdr_cutoff"],
        lfc                  = config["de"]["log2fc_cutoff"],
        tumor_label          = config["de"]["tumor_label"],
        normal_label         = config["de"]["normal_label"],
        de_sig_by            = config["de"].get("de_sig_by", "pvalue"),
        heatmap_top_n        = config["de"].get("heatmap_top_n", 10),
        fsj_concordance_lfc  = config["de"].get("fsj_concordance_lfc", 0.0),
        filter_expr_min_count       = config["de"].get("filter_expr_min_count", 5),
        filter_expr_min_total_count = config["de"].get("filter_expr_min_total_count", 15),
        deseq2_size_factor_type     = config["de"].get("deseq2_size_factor_type", "poscounts"),
        paired_design_override      = config["de"].get("paired_design_override", "auto"),
    log: "logs/de_analysis.log"
    script:
        "../../scripts/analysis.R"


rule annotate_circbase:
    """Annotate all consensus circRNAs (merged across samples) against circBase hg19."""
    input:
        summary = expand(
            RESULTS_DIR + "/circRNA/{srr}/consensus_summary.tsv", srr=SAMPLES
        ),
    output:
        RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
    params:
        slop   = config["consensus"]["slop"],
        cb_arg = lambda w: (
            f"--circbase-file {config['circbase_file']}"
            if config.get("circbase_file") else ""
        ),
    log: "logs/circbase_annotation.log"
    shell:
        """
        python scripts/annotate_circbase.py \
            --summary     {input.summary} \
            --output      {output} \
            --slop        {params.slop} \
            {params.cb_arg} \
            2>&1 | tee {log}
        """


rule rank_biomarkers:
    """Rank DE circRNAs by composite biomarker score (6D when interactions available) + n_sig_methods."""
    input:
        de           = RESULTS_DIR + "/de/de_results.tsv",
        de_deseq2    = RESULTS_DIR + "/de/de_results_deseq2.tsv",
        de_limma     = RESULTS_DIR + "/de/de_results_limma.tsv",
        annot        = RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
        summary      = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        interactions = RESULTS_DIR + "/de/interactions.json",
    output:
        RESULTS_DIR + "/de/biomarker_candidates.tsv",
        # NOTE: biomarker_candidates_fully_evaluated.tsv (--fully-evaluated-only)
        # is intentionally NOT declared as a rule output yet. Declaring a new
        # output here would change this rule's output signature and make
        # Snakemake consider rank_biomarkers "incomplete" for every already-
        # completed dataset the next time its DAG is touched, forcing a rerun
        # of this step across all 19 published datasets -- which contradicts
        # the explicit "don't touch existing results yet" decision (2026-10).
        # Run the script directly (see CLAUDE.md) when the fully-evaluated
        # table is needed for a specific dataset; wire it into this rule only
        # once a deliberate, dataset-by-dataset rerun has been decided.
    params:
        fdr             = config["de"]["fdr_cutoff"],
        lfc             = config["de"]["log2fc_cutoff"],
        de_sig_by       = config["de"].get("de_sig_by", "pvalue"),
        # Only emit --sig-cap/--fc-cap when they differ from rank_biomarkers.py's
        # own defaults (10 / 5), to keep the default-config command byte-identical.
        cap_args = " ".join(filter(None, [
            f"--sig-cap {config['de'].get('biomarker', {}).get('sig_cap', 10)}"
            if config["de"].get("biomarker", {}).get("sig_cap", 10) != 10 else "",
            f"--fc-cap {config['de'].get('biomarker', {}).get('fc_cap', 5)}"
            if config["de"].get("biomarker", {}).get("fc_cap", 5) != 5 else "",
        ])),
        # Only emit --components-json when weights differ from the equal-weight
        # default, so the default-config command stays byte-identical.
        components_arg = (
            lambda w: (
                "--components-json " + __import__("shlex").quote(__import__("json").dumps(
                    config["de"]["biomarker"]["components"]))
                if config["de"].get("biomarker", {}).get("components")
                   != {k: {"enabled": True, "weight": 1.0} for k in
                       ("sig", "fc", "conf", "known", "mirna", "rbp")}
                else ""
            )
        ),
    log: "logs/rank_biomarkers.log"
    shell:
        """
        python scripts/rank_biomarkers.py \
            --de           {input.de} \
            --de-deseq2    {input.de_deseq2} \
            --de-limma     {input.de_limma} \
            --annot        {input.annot} \
            --output       {output} \
            --interactions {input.interactions} \
            --fdr     {params.fdr} \
            --lfc     {params.lfc} \
            --de-sig-by {params.de_sig_by} \
            {params.cap_args} \
            {params.components_arg} \
            2>&1 | tee {log}
        """


rule assign_isoforms:
    """Map consensus circRNA coordinates to host genes; build isoform group table."""
    input:
        matrix = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        gtf    = config["genome"]["gtf"],
    output:
        RESULTS_DIR + "/circRNA/isoform_groups.tsv",
    log: "logs/assign_isoforms.log"
    shell:
        """
        python scripts/assign_isoforms.py \
            --count-matrix {input.matrix} \
            --gtf          {input.gtf} \
            --out          {output} \
            2>&1 | tee {log}
        """


rule isoform_switching:
    """Compute per-isoform IUI and test for switching between tumor and normal."""
    input:
        bsj_matrix     = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        isoform_groups = RESULTS_DIR + "/circRNA/isoform_groups.tsv",
        de_results     = RESULTS_DIR + "/de/de_results.tsv",
        sample_groups  = config["groups"],
    output:
        iui_matrix = RESULTS_DIR + "/de/iui_matrix.tsv",
        switching  = RESULTS_DIR + "/de/isoform_switching.tsv",
    params:
        tumor_label      = config["de"]["tumor_label"],
        normal_label     = config["de"]["normal_label"],
        fdr              = config["de"].get("isoform_fdr_cutoff", 0.1),
        delta_iui_cutoff = config["de"].get("delta_iui_cutoff", 0.1),
    log: "logs/isoform_switching.log"
    script:
        "../../scripts/isoform_switching.R"


rule predict_interactions:
    """Query CircInteractome for miRNA/RBP binding sites of top DE circRNAs (union of all 3 methods)."""
    input:
        de       = RESULTS_DIR + "/de/de_results.tsv",
        de_edger = RESULTS_DIR + "/de/de_results_edgeR_ciriquant.tsv",
        de_deseq = RESULTS_DIR + "/de/de_results_deseq2.tsv",
        de_limma = RESULTS_DIR + "/de/de_results_limma.tsv",
        iso      = RESULTS_DIR + "/circRNA/isoform_groups.tsv",
        circbase = RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
    output:
        RESULTS_DIR + "/de/interactions.json",
    params:
        top_n        = config["de"].get("interaction_top_n", 50),
        clip_exp_num = config["de"].get("encori_clip_exp_num", 1),
        program_num  = config["de"].get("encori_program_num", 2),
        gtf          = config["genome"]["gtf"],
    log: "logs/predict_interactions.log"
    shell:
        """
        python scripts/predict_interactions.py \
            --de           {input.de} \
            --iso          {input.iso} \
            --circbase     {input.circbase} \
            --out          {output} \
            --gtf          {params.gtf} \
            --de-edger     {input.de_edger} \
            --de-deseq2    {input.de_deseq} \
            --de-limma     {input.de_limma} \
            --top-n        {params.top_n} \
            --clip-exp-num {params.clip_exp_num} \
            --program-num  {params.program_num} \
            2>&1 | tee {log}
        """


# ─── Effect-size stability (leave-one-sample-out) and per-method candidate tables ───
# Added 2026-10: the LOO fits are dataset-wide and shared by all three primary methods.
# New outputs only (no change to the rank_biomarkers rule above), so published datasets
# are not scheduled for upstream reruns; see the dry-run record in audit/.
import csv as _csv
_LOO_SAMPLES = [r["srr_id"] for r in _csv.DictReader(open(config["groups"], encoding="utf-8"))]
_LOO_OUT = ([f"{RESULTS_DIR}/de/loo/{d}/{m}.tsv" for d in ["NONE"] + _LOO_SAMPLES
             for m in ("edger", "deseq2", "limma")]
            + [f"{RESULTS_DIR}/de/loo/run_info.tsv"])
_BM_V2_PRIMARY = {"edgeR": "de_results.tsv", "DESeq2": "de_results_deseq2.tsv",
                  "limma": "de_results_limma.tsv"}


rule loo_stability:
    """Leave-one-sample-out runs of the three DE chains (mirrors analysis.R; paired if patient_id)."""
    input:
        matrix = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        fsj    = RESULTS_DIR + "/circRNA/fsj_count_matrix.tsv",
        groups = config["groups"],
    output:
        _LOO_OUT
    params:
        results_dir  = RESULTS_DIR,
        out_dir      = RESULTS_DIR + "/de/loo",
        tumor        = config["de"]["tumor_label"],
        normal       = config["de"]["normal_label"],
    log: "logs/loo_stability.log"
    shell:
        "Rscript scripts/loo_stability.R {params.results_dir} {input.groups} {params.out_dir} "
        "{params.tumor} {params.normal} > {log} 2>&1"


rule rank_biomarkers_v2:
    """Candidate table for one primary method: cross-method flag + LOO stability; default view = passing."""
    input:
        de_edger     = RESULTS_DIR + "/de/de_results.tsv",
        de_deseq2    = RESULTS_DIR + "/de/de_results_deseq2.tsv",
        de_limma     = RESULTS_DIR + "/de/de_results_limma.tsv",
        annot        = RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
        isoforms     = RESULTS_DIR + "/circRNA/isoform_groups.tsv",
        interactions = RESULTS_DIR + "/de/interactions.json",
        loo          = _LOO_OUT,
    output:
        candidates = RESULTS_DIR + "/de/biomarker_candidates_v2_{method}.tsv",
        status     = RESULTS_DIR + "/de/biomarker_candidates_v2_{method}.status.json",
    wildcard_constraints:
        method = "edgeR|DESeq2|limma"
    params:
        results_dir         = RESULTS_DIR,
        primary             = lambda w: _BM_V2_PRIMARY[w.method],
        fdr                 = config["de"]["fdr_cutoff"],
        lfc                 = config["de"]["log2fc_cutoff"],
        de_sig_by           = config["de"].get("de_sig_by", "pvalue"),
        cross_method_min    = config["de"].get("biomarker", {}).get("cross_method_min", 3),
        stability_min_ratio = config["de"].get("biomarker", {}).get("stability_min_ratio", 0.5),
        loo_dir             = RESULTS_DIR + "/de/loo",
        # Same wiring as the rank_biomarkers rule: emit --components-json only when the
        # weights differ from the equal-weight default (default output is unchanged).
        components_arg = (
            lambda w: (
                "--components-json " + __import__("shlex").quote(__import__("json").dumps(
                    config["de"]["biomarker"]["components"]))
                if config["de"].get("biomarker", {}).get("components")
                   != {k: {"enabled": True, "weight": 1.0} for k in
                       ("sig", "fc", "conf", "known", "mirna", "rbp")}
                else ""
            )
        ),
        # Same wiring as the v1 rank_biomarkers rule: emit each flag only when it differs from
        # rank_biomarkers.py's default, so the default command is unchanged.
        sig_cap_arg = (f"--sig-cap {config['de'].get('biomarker', {}).get('sig_cap', 10)}"
                       if config["de"].get("biomarker", {}).get("sig_cap", 10) != 10 else ""),
        fc_cap_arg  = (f"--fc-cap {config['de'].get('biomarker', {}).get('fc_cap', 5)}"
                       if config["de"].get("biomarker", {}).get("fc_cap", 5) != 5 else ""),
        stability_max_arg = (f"--stability-max-candidates "
                             f"{config['de'].get('biomarker', {}).get('stability_max_candidates', 500)}"
                             if config["de"].get("biomarker", {}).get("stability_max_candidates", 500) != 500
                             else ""),
    log: "logs/rank_biomarkers_v2_{method}.log"
    shell:
        "python scripts/rank_biomarkers.py --de {params.results_dir}/de/{params.primary} "
        "--primary-method {wildcards.method} --de-edger {input.de_edger} "
        "--de-deseq2 {input.de_deseq2} --de-limma {input.de_limma} "
        "--annot {input.annot} --isoforms {input.isoforms} --interactions {input.interactions} "
        "--output {output.candidates} --status-json {output.status} "
        "--fdr {params.fdr} --lfc {params.lfc} --de-sig-by {params.de_sig_by} "
        "--loo-dir {params.loo_dir} --cross-method-min {params.cross_method_min} "
        "--stability-min-ratio {params.stability_min_ratio} {params.components_arg} "
        "{params.sig_cap_arg} {params.fc_cap_arg} {params.stability_max_arg} > {log} 2>&1"


rule generate_report:
    """Build a self-contained HTML summary report."""
    input:
        de       = RESULTS_DIR + "/de/de_results.tsv",
        de_edger = RESULTS_DIR + "/de/de_results_edgeR_ciriquant.tsv",
        de_deseq = RESULTS_DIR + "/de/de_results_deseq2.tsv",
        de_limma = RESULTS_DIR + "/de/de_results_limma.tsv",
        biomarkers_edger  = RESULTS_DIR + "/de/biomarker_candidates_v2_edgeR.tsv",
        biomarkers_deseq2 = RESULTS_DIR + "/de/biomarker_candidates_v2_DESeq2.tsv",
        biomarkers_limma  = RESULTS_DIR + "/de/biomarker_candidates_v2_limma.tsv",
        biomarker_status  = expand(RESULTS_DIR + "/de/biomarker_candidates_v2_{m}.status.json",
                                   m=["edgeR", "DESeq2", "limma"]),
        matrix         = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        volcano        = RESULTS_DIR + "/plots/volcano.pdf",
        heatmap        = RESULTS_DIR + "/plots/heatmap.pdf",
        pca            = RESULTS_DIR + "/plots/pca.pdf",
        multiqc        = RESULTS_DIR + "/qc/multiqc_report.html",
        switching      = RESULTS_DIR + "/de/isoform_switching.tsv",
        groups         = config["groups"],
        isoform_groups = RESULTS_DIR + "/circRNA/isoform_groups.tsv",
        circbase_annot = RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
        interactions   = RESULTS_DIR + "/de/interactions.json",
    output:
        RESULTS_DIR + "/report.html",
    params:
        project_id   = config["project_id"],
        fdr          = config["de"]["fdr_cutoff"],
        lfc          = config["de"]["log2fc_cutoff"],
        de_method    = DE_METHOD,
        de_sig_by     = config["de"].get("de_sig_by", "pvalue"),
        tumor_label   = config["de"]["tumor_label"],
        normal_label  = config["de"]["normal_label"],
        heatmap_top_n = config["de"].get("heatmap_top_n", 10),
        study_title   = config.get("study_title", ""),
    log: "logs/generate_report.log"
    script:
        "../../scripts/generate_report.py"
