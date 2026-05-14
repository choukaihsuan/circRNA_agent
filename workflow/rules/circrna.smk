"""
circRNA detection rules:
  fastp-trimmed FASTQ → CIRIquant + STAR/DCC → consensus filter → count matrix
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
        r1  = TRIMMED_DIR + "/{srr}_1.fastq.gz",
        r2  = TRIMMED_DIR + "/{srr}_2.fastq.gz",
        cfg = config["ciriquant_config"],
        _   = "config/.ciriquant_ready",
    output:
        gtf = "results/circRNA/{srr}/{srr}.gtf",
        bsj = "results/circRNA/{srr}/{srr}.bsj",
    params:
        outdir = "results/circRNA/{srr}",
        sample = "{srr}",
    threads: config["threads"]
    resources:
        mem_gb = 16
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


rule star_align:
    """STAR alignment for DCC chimeric junction detection."""
    input:
        r1 = TRIMMED_DIR + "/{srr}_1.fastq.gz",
        r2 = TRIMMED_DIR + "/{srr}_2.fastq.gz",
    output:
        junc = "results/circRNA/{srr}/Chimeric.out.junction",
        bam  = "results/circRNA/{srr}/Aligned.sortedByCoord.out.bam",
    params:
        index   = config["genome"]["star_index"],
        prefix  = "results/circRNA/{srr}/",
        tmp_dir = "/home/choukaihsuan/star_tmp/{srr}",
    threads: config["threads"]
    resources:
        mem_gb = 26
    log: "logs/star/{srr}.log"
    shell:
        """
        mkdir -p {params.tmp_dir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmp_dir}/_STARtmp \
            --chimSegmentMin 10 \
            --chimOutType Junctions \
            --alignSJDBoverhangMin 10 \
            --outSAMattributes NH HI NM MD AS \
            > {log} 2>&1
        samtools index {output.bam}
        rm -rf {params.tmp_dir}
        """


rule dcc:
    """Detect circRNAs with DCC using STAR chimeric junctions."""
    input:
        junc = "results/circRNA/{srr}/Chimeric.out.junction",
        bam  = "results/circRNA/{srr}/Aligned.sortedByCoord.out.bam",
        gtf  = config["genome"]["gtf"],
    output:
        "results/circRNA/{srr}/DCC/CircCoordinates",
    params:
        outdir = "results/circRNA/{srr}/DCC",
    threads: 4
    log: "logs/dcc/{srr}.log"
    shell:
        """
        mkdir -p {params.outdir}
        echo "{input.junc}" > {params.outdir}/junction.list
        dcc @{params.outdir}/junction.list \
            -mt1 @{params.outdir}/junction.list \
            -D -an {input.gtf} \
            -Pi -F -M -Nr 5 1 \
            -fg -G -A {input.bam} \
            -O {params.outdir} \
            -T {threads} \
            > {log} 2>&1
        """


rule consensus_filter:
    """Vote across CIRIquant and DCC; keep circRNAs supported by both tools."""
    input:
        cirique = "results/circRNA/{srr}/{srr}.gtf",
        dcc     = "results/circRNA/{srr}/DCC/CircCoordinates",
    output:
        bed     = "results/circRNA/{srr}/high_confidence.bed",
        summary = "results/circRNA/{srr}/consensus_summary.tsv",
    params:
        min_tools = config["consensus"]["min_tools"],
        slop      = config["consensus"]["slop"],
        min_bsj   = config["consensus"]["min_bsj_reads"],
    log: "logs/consensus/{srr}.log"
    shell:
        """
        python scripts/consensus_filter.py \
            --cirique   {input.cirique} \
            --dcc       {input.dcc} \
            --output    {output.bed} \
            --summary   {output.summary} \
            --min-tools {params.min_tools} \
            --slop      {params.slop} \
            --min-bsj   {params.min_bsj} \
            2>&1 | tee {log}
        """


rule merge_counts:
    """Filter by consensus BED and build circRNA × sample BSJ count matrix."""
    input:
        gtfs = expand("results/circRNA/{srr}/{srr}.gtf", srr=SAMPLES),
        beds = expand("results/circRNA/{srr}/high_confidence.bed", srr=SAMPLES),
    output:
        matrix = "results/circRNA/count_matrix.tsv",
    shell:
        """
        python scripts/merge_counts.py \
            --gtfs       {input.gtfs} \
            --output     {output.matrix} \
            --filter-bed {input.beds}
        """
