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


rule generate_report:
    """Build a self-contained HTML summary report."""
    input:
        de      = RESULTS_DIR + "/de/de_results.tsv",
        matrix  = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        volcano = RESULTS_DIR + "/plots/volcano.pdf",
        heatmap = RESULTS_DIR + "/plots/heatmap.pdf",
        pca     = RESULTS_DIR + "/plots/pca.pdf",
        multiqc = RESULTS_DIR + "/qc/multiqc_report.html",
    output:
        RESULTS_DIR + "/report.html",
    params:
        project_id = config["project_id"],
        fdr        = config["de"]["fdr_cutoff"],
        lfc        = config["de"]["log2fc_cutoff"],
    log: "logs/generate_report.log"
    script:
        "../../scripts/generate_report.py"
