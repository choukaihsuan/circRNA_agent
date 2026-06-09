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
        de_sig_by            = config["de"].get("de_sig_by", "auto"),
        heatmap_top_n        = config["de"].get("heatmap_top_n", 10),
        fsj_concordance_lfc  = config["de"].get("fsj_concordance_lfc", 0.0),
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
    params:
        fdr             = config["de"]["fdr_cutoff"],
        lfc             = config["de"]["log2fc_cutoff"],
        de_sig_by       = config["de"].get("de_sig_by", "auto"),
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
    """Query CircInteractome for miRNA/RBP binding sites of top DE circRNAs."""
    input:
        de       = RESULTS_DIR + "/de/de_results.tsv",
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
            --top-n        {params.top_n} \
            --clip-exp-num {params.clip_exp_num} \
            --program-num  {params.program_num} \
            2>&1 | tee {log}
        """


rule generate_report:
    """Build a self-contained HTML summary report."""
    input:
        de       = RESULTS_DIR + "/de/de_results.tsv",
        de_edger = RESULTS_DIR + "/de/de_results_edgeR_ciriquant.tsv",
        de_deseq = RESULTS_DIR + "/de/de_results_deseq2.tsv",
        de_limma = RESULTS_DIR + "/de/de_results_limma.tsv",
        biomarkers     = RESULTS_DIR + "/de/biomarker_candidates.tsv",
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
        de_sig_by     = config["de"].get("de_sig_by", "auto"),
        tumor_label   = config["de"]["tumor_label"],
        normal_label  = config["de"]["normal_label"],
        heatmap_top_n = config["de"].get("heatmap_top_n", 10),
    log: "logs/generate_report.log"
    script:
        "../../scripts/generate_report.py"
