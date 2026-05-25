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
    """Annotate consensus circRNAs against the circBase hg19 database."""
    input:
        summary = expand(
            RESULTS_DIR + "/circRNA/{srr}/consensus_summary.tsv", srr=SAMPLES
        ),
    output:
        RESULTS_DIR + "/circRNA/circbase_annotated.tsv",
    params:
        slop           = config["consensus"]["slop"],
        circbase_file  = config.get("circbase_file", ""),
    log: "logs/circbase_annotation.log"
    run:
        import subprocess, sys
        summaries = " ".join(input.summary)
        cb_arg = f"--circbase-file {params.circbase_file}" if params.circbase_file else ""
        # Merge all per-sample summaries first, then annotate the first one
        # (annotation is per-coordinate, sample-agnostic)
        # We annotate the merged count matrix coordinates via the first summary
        shell(
            f"python scripts/annotate_circbase.py "
            f"--summary {input.summary[0]} "
            f"--output {output[0]} "
            f"--slop {params.slop} "
            f"{cb_arg} "
            f"2>&1 | tee {log}"
        )


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
