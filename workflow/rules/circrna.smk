"""
circRNA detection rules:
  fastp-trimmed FASTQ → CIRIquant + STAR/DCC → consensus filter → count matrix
"""

import os


wildcard_constraints:
    srr = r"[A-Z]+\d+"


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
    """Run CIRIquant on one sample; outputs per-sample GTF."""
    input:
        r1  = TRIMMED_DIR + "/{srr}_1.fastq.gz",
        r2  = TRIMMED_DIR + "/{srr}_2.fastq.gz",
        cfg = config["ciriquant_config"],
        _   = "config/.ciriquant_ready",
    output:
        gtf = RESULTS_DIR + "/circRNA/{srr}/{srr}.gtf",
    params:
        outdir = RESULTS_DIR + "/circRNA/{srr}",
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
        rm -rf {params.outdir}/circ {params.outdir}/align {params.outdir}/gene
        """


rule star_align:
    """STAR paired-end alignment for DCC chimeric junction detection."""
    input:
        r1 = TRIMMED_DIR + "/{srr}_1.fastq.gz",
        r2 = TRIMMED_DIR + "/{srr}_2.fastq.gz",
    output:
        junc = RESULTS_DIR + "/circRNA/{srr}/Chimeric.out.junction",
        bam  = RESULTS_DIR + "/circRNA/{srr}/Aligned.sortedByCoord.out.bam",
    params:
        index   = config["genome"]["star_index"],
        prefix  = RESULTS_DIR + "/circRNA/{srr}/",
        tmp_dir = RESULTS_DIR + "/circRNA/{srr}/star_tmp",
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


rule star_align_mate1:
    """STAR single-end alignment of mate 1 only (required by DCC -mt1)."""
    input:
        r1 = TRIMMED_DIR + "/{srr}_1.fastq.gz",
    output:
        junc = RESULTS_DIR + "/circRNA/{srr}/mate1/Chimeric.out.junction",
    params:
        index   = config["genome"]["star_index"],
        prefix  = RESULTS_DIR + "/circRNA/{srr}/mate1/",
        tmp_dir = RESULTS_DIR + "/circRNA/{srr}/mate1_tmp",
    threads: config["threads"]
    resources:
        mem_gb = 26
    log: "logs/star/{srr}_mate1.log"
    shell:
        """
        mkdir -p {params.tmp_dir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.r1} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmp_dir}/_STARtmp \
            --chimSegmentMin 10 \
            --chimOutType Junctions \
            --alignSJDBoverhangMin 10 \
            > {log} 2>&1
        rm -rf {params.tmp_dir}
        """


rule star_align_mate2:
    """STAR single-end alignment of mate 2 only (required by DCC -mt2)."""
    input:
        r2 = TRIMMED_DIR + "/{srr}_2.fastq.gz",
    output:
        junc = RESULTS_DIR + "/circRNA/{srr}/mate2/Chimeric.out.junction",
    params:
        index   = config["genome"]["star_index"],
        prefix  = RESULTS_DIR + "/circRNA/{srr}/mate2/",
        tmp_dir = RESULTS_DIR + "/circRNA/{srr}/mate2_tmp",
    threads: config["threads"]
    resources:
        mem_gb = 26
    log: "logs/star/{srr}_mate2.log"
    shell:
        """
        mkdir -p {params.tmp_dir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmp_dir}/_STARtmp \
            --chimSegmentMin 10 \
            --chimOutType Junctions \
            --alignSJDBoverhangMin 10 \
            > {log} 2>&1
        rm -rf {params.tmp_dir}
        """


rule dcc:
    """Detect circRNAs with DCC using STAR chimeric junctions (paired + mate1 + mate2)."""
    input:
        junc  = RESULTS_DIR + "/circRNA/{srr}/Chimeric.out.junction",
        junc1 = RESULTS_DIR + "/circRNA/{srr}/mate1/Chimeric.out.junction",
        junc2 = RESULTS_DIR + "/circRNA/{srr}/mate2/Chimeric.out.junction",
        bam   = RESULTS_DIR + "/circRNA/{srr}/Aligned.sortedByCoord.out.bam",
        gtf   = config["genome"]["gtf"],
    output:
        RESULTS_DIR + "/circRNA/{srr}/DCC/CircCoordinates",
    params:
        outdir = RESULTS_DIR + "/circRNA/{srr}/DCC",
    threads: 4
    log: "logs/dcc/{srr}.log"
    shell:
        """
        mkdir -p {params.outdir}
        ( cd {params.outdir} && DCC {input.junc} \
            -mt1 {input.junc1} \
            -mt2 {input.junc2} \
            -D -an {input.gtf} \
            -Pi -F -M -Nr 5 1 \
            -G -B {input.bam} \
            -O {params.outdir} \
            -T {threads} ) \
            > {log} 2>&1
        """


rule consensus_filter:
    """Filter circRNAs; supports CIRIquant-only, DCC-only, or consensus mode."""
    input:
        cirique = (RESULTS_DIR + "/circRNA/{srr}/{srr}.gtf")
                  if USE_CIRIQUANT else [],
        dcc     = (RESULTS_DIR + "/circRNA/{srr}/DCC/CircCoordinates")
                  if USE_DCC else [],
    output:
        bed     = RESULTS_DIR + "/circRNA/{srr}/high_confidence.bed",
        summary = RESULTS_DIR + "/circRNA/{srr}/consensus_summary.tsv",
    params:
        cirique_arg       = lambda w, input: f"--cirique {input.cirique}" if USE_CIRIQUANT else "",
        dcc_arg           = lambda w, input: f"--dcc {input.dcc}"         if USE_DCC       else "",
        min_tools         = config["consensus"]["min_tools"],
        slop              = config["consensus"]["slop"],
        min_bsj           = config["consensus"]["min_bsj_reads"],
        max_junction_ratio = config["consensus"].get("max_junction_ratio", 1.0),
    log: "logs/consensus/{srr}.log"
    shell:
        """
        python scripts/consensus_filter.py \
            {params.cirique_arg} \
            {params.dcc_arg} \
            --output              {output.bed} \
            --summary             {output.summary} \
            --min-tools           {params.min_tools} \
            --slop                {params.slop} \
            --min-bsj             {params.min_bsj} \
            --max-junction-ratio  {params.max_junction_ratio} \
            2>&1 | tee {log}
        """


rule merge_counts:
    """Filter by consensus BED and build circRNA × sample BSJ + FSJ count matrices."""
    input:
        gtfs = expand(RESULTS_DIR + "/circRNA/{srr}/{srr}.gtf", srr=SAMPLES),
        beds = expand(RESULTS_DIR + "/circRNA/{srr}/high_confidence.bed", srr=SAMPLES),
    output:
        matrix     = RESULTS_DIR + "/circRNA/count_matrix.tsv",
        fsj_matrix = RESULTS_DIR + "/circRNA/fsj_count_matrix.tsv",
    shell:
        """
        python scripts/merge_counts.py \
            --gtfs       {input.gtfs} \
            --output     {output.matrix} \
            --output-fsj {output.fsj_matrix} \
            --filter-bed {input.beds}
        """
