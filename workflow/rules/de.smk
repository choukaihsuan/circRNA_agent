"""
Differential expression and reporting rules (requires sample_groups.csv).
Skipped automatically when HAS_GROUPS is False (see Snakefile).
"""


rule de_analysis:
    """DESeq2-based DE analysis of circRNAs between tumor and normal groups."""
    input:
        matrix = "results/circRNA/count_matrix.tsv",
        groups = config["groups"],
    output:
        de      = "results/de/de_results.tsv",
        volcano = "results/plots/volcano.pdf",
        heatmap = "results/plots/heatmap.pdf",
        pca     = "results/plots/pca.pdf",
    params:
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
        de      = "results/de/de_results.tsv",
        matrix  = "results/circRNA/count_matrix.tsv",
        volcano = "results/plots/volcano.pdf",
        heatmap = "results/plots/heatmap.pdf",
        pca     = "results/plots/pca.pdf",
        multiqc = "results/qc/multiqc_report.html",
    output:
        "results/report.html",
    params:
        project_id = config["project_id"],
    script:
        "../../scripts/generate_report.py"
