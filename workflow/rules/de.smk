"""
Differential expression and reporting rules (requires sample_groups.csv).
Skipped automatically when HAS_GROUPS is False (see Snakefile).
"""


rule de_analysis:
    """circRNA DE analysis; method = edgeR_ciriquant / deseq2 / limma (config de.method)."""
    input:
        matrix     = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        fsj_matrix = RESULTS_DIR + "/circRNA/fsj_count_matrix.tsv",
        groups     = config["groups"],
    output:
        de      = RESULTS_DIR + "/de/de_results.tsv",
        volcano = RESULTS_DIR + "/plots/volcano.pdf",
        heatmap = RESULTS_DIR + "/plots/heatmap.pdf",
        pca     = RESULTS_DIR + "/plots/pca.pdf",
    params:
        de_method    = DE_METHOD,
        fdr          = config["de"]["fdr_cutoff"],
        lfc          = config["de"]["log2fc_cutoff"],
        tumor_label  = config["de"]["tumor_label"],
        normal_label = config["de"]["normal_label"],
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
    """Rank DE circRNAs by composite biomarker score."""
    input:
        de      = RESULTS_DIR + "/de/de_results.tsv",
        annot   = RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
        summary = RESULTS_DIR + "/circRNA/count_matrix.tsv",
    output:
        RESULTS_DIR + "/de/biomarker_candidates.tsv",
    params:
        fdr = config["de"]["fdr_cutoff"],
        lfc = config["de"]["log2fc_cutoff"],
    log: "logs/rank_biomarkers.log"
    shell:
        """
        python scripts/rank_biomarkers.py \
            --de      {input.de} \
            --annot   {input.annot} \
            --output  {output} \
            --fdr     {params.fdr} \
            --lfc     {params.lfc} \
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


rule generate_report:
    """Build a self-contained HTML summary report."""
    input:
        de         = RESULTS_DIR + "/de/de_results.tsv",
        biomarkers = RESULTS_DIR + "/de/biomarker_candidates.tsv",
        matrix     = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        volcano    = RESULTS_DIR + "/plots/volcano.pdf",
        heatmap    = RESULTS_DIR + "/plots/heatmap.pdf",
        pca        = RESULTS_DIR + "/plots/pca.pdf",
        multiqc    = RESULTS_DIR + "/qc/multiqc_report.html",
        switching  = RESULTS_DIR + "/de/isoform_switching.tsv",
        groups     = config["groups"],
    output:
        RESULTS_DIR + "/report.html",
    params:
        project_id = config["project_id"],
        fdr        = config["de"]["fdr_cutoff"],
        lfc        = config["de"]["log2fc_cutoff"],
        de_method  = DE_METHOD,
    log: "logs/generate_report.log"
    script:
        "../../scripts/generate_report.py"
