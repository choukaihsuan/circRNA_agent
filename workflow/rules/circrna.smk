"""
circRNA detection rules:
  fastp-trimmed FASTQ → CIRIquant (per sample) → merged count matrix
"""

import os


rule check_ciriquant_config:
    """Validate that the CIRIquant config YAML exists before any sample runs."""
    output:
        touch("config/.ciriquant_ready"),
    run:
        cfg = config["ciriquant_config"]
        if not os.path.exists(cfg):
            raise ValueError(
                f"CIRIquant config not found: {cfg}\n"
                "Generate it with:  python scripts/agent.py --setup-ciriquant"
            )


rule ciriquant:
    """Run CIRIquant on one sample; outputs per-sample GTF and BSJ table."""
    input:
        r1  = "data/trimmed/{srr}_1.fastq.gz",
        r2  = "data/trimmed/{srr}_2.fastq.gz",
        cfg = config["ciriquant_config"],
        _   = "config/.ciriquant_ready",
    output:
        gtf = "results/circRNA/{srr}/{srr}.gtf",
        bsj = "results/circRNA/{srr}/{srr}.bsj",
    params:
        outdir = "results/circRNA/{srr}",
        sample = "{srr}",
    threads: config["threads"]
    log: "logs/ciriquant/{srr}.log"
    shell:
        """
        CIRIquant \
            -t {threads} \
            -1 {input.r1} \
            -2 {input.r2} \
            --config {input.cfg} \
            -o {params.outdir} \
            -p {params.sample} \
            > {log} 2>&1
        """


rule merge_counts:
    """Parse all per-sample GTF files and build a sample × circRNA BSJ count matrix."""
    input:
        gtfs = expand("results/circRNA/{srr}/{srr}.gtf", srr=SAMPLES),
    output:
        matrix = "results/circRNA/count_matrix.tsv",
    script:
        "../../scripts/merge_counts.py"
