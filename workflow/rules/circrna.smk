"""
circRNA detection rules:
  fastp-trimmed FASTQ → CIRIquant + STAR/DCC → consensus filter → count matrix
"""

import os


# disk_gb resource claims below (ciriquant, star_align*) are a single shared
# ESTIMATE (200 GB/sample), not a measured value. Provenance: a prior
# benchmark run (SRR444655, large single total-RNA sample) filled the disk
# mid-run (`OSError: [Errno 28] No space left on device`); the real
# per-sample peak was never captured by a df/du time series. 200 GB is
# reverse-engineered from two indirect signals: (a) the run completed
# successfully once ~600GB of unrelated data was freed, and (b) this file's
# own `ciriquant` rule already `rm -rf`s its circ/align/gene subdirectories
# after the tool finishes, implying the original author already knew those
# intermediates were large enough to be worth deleting. No rule here differs
# in intermediate-file strategy between samples/datasets, so the same number
# is applied to ciriquant and all three star_align* rules pending an actual
# measurement -- see logs/disk_usage_{srr}.tsv (written by this file's
# disk-sampler, scripts/sample_disk_usage.py) for real data once a future
# run's samples accumulate it.
_DISK_GB_PER_SAMPLE_ESTIMATE = 200

wildcard_constraints:
    srr = r"[A-Z]+\d+"


rule check_ciriquant_config:
    """Validate that the CIRIquant config YAML exists before any sample runs."""
    output:
        touch("config/.ciriquant_ready"),
    run:
        cfg = CIRIQUANT_CONFIG
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
        cfg = CIRIQUANT_CONFIG,
        _   = "config/.ciriquant_ready",
    output:
        gtf = RESULTS_DIR + "/circRNA/{srr}/{srr}.gtf",
    params:
        outdir   = RESULTS_DIR + "/circRNA/{srr}",
        sample   = "{srr}",
        # Only emit -a when it differs from CIRIquant's own default (5), so the
        # generated command at default settings is byte-identical to the
        # previously-hardcoded version (see tests/test_validate_config.py).
        anchor_arg = (f"-a {config['tools']['ciriquant']['anchor']}"
                      if config["tools"]["ciriquant"]["anchor"] != 5 else ""),
        no_gene  = "--no-gene" if config["tools"]["ciriquant"]["no_gene"] else "",
    threads: config["threads"]
    resources:
        mem_gb = 16,
        disk_gb = _DISK_GB_PER_SAMPLE_ESTIMATE
    log: "logs/ciriquant/{srr}.log"
    # Background disk sampler below (mkdir/python.../kill lines): near-zero
    # cost (one df + one `du -sb` per minute), piggybacks on whatever real
    # sample happens to run next through this rule -- never run this rule, or
    # a dedicated analysis, just to collect this data (see
    # scripts/sample_disk_usage.py's own docstring). Writes
    # logs/disk_usage_{srr}.tsv; killed unconditionally after CIRIquant
    # finishes (success or failure) so it never outlives the rule.
    shell:
        """
        mkdir -p {params.outdir}
        python scripts/sample_disk_usage.py {params.outdir} logs/disk_usage_{wildcards.srr}.tsv &
        _SAMPLER_PID=$!
        CIRIquant \
            -t {threads} \
            -1 {input.r1} \
            -2 {input.r2} \
            --config {input.cfg} \
            -o {params.outdir} \
            -p {params.sample} \
            {params.anchor_arg} \
            {params.no_gene} \
            > {log} 2>&1
        kill $_SAMPLER_PID 2>/dev/null || true
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
        mem_gb = 26,
        disk_gb = _DISK_GB_PER_SAMPLE_ESTIMATE
    log: "logs/star/{srr}.log"
    shell:
        """
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM Unsorted \
            --outFileNamePrefix {params.prefix} \
            --outTmpDir {params.tmp_dir}/_STARtmp \
            --chimSegmentMin 10 \
            --chimOutType Junctions \
            --alignSJDBoverhangMin 10 \
            --outSAMattributes NH HI NM MD AS \
            > {log} 2>&1
        samtools sort -@ {threads} -m 2G \
            -o {output.bam} \
            {params.prefix}Aligned.out.bam \
            >> {log} 2>&1
        rm -f {params.prefix}Aligned.out.bam
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
        mem_gb = 26,
        disk_gb = _DISK_GB_PER_SAMPLE_ESTIMATE
    log: "logs/star/{srr}_mate1.log"
    shell:
        """
        ulimit -n 65536 || true
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
        mem_gb = 26,
        disk_gb = _DISK_GB_PER_SAMPLE_ESTIMATE
    log: "logs/star/{srr}_mate2.log"
    shell:
        """
        ulimit -n 65536 || true
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
        outdir   = RESULTS_DIR + "/circRNA/{srr}/DCC",
        nr       = lambda w: f"{config['tools']['dcc']['count_threshold']} {config['tools']['dcc']['replicate_threshold']}",
        filter_f = "-F" if config["tools"]["dcc"]["filter_recommended"] else "",
        filter_m = "-M" if config["tools"]["dcc"]["filter_mito"] else "",
        # Only emit -E/-m/-n when they differ from DCC's own defaults
        # (5 / 1,000,000 / 30), so the default-config command is byte-identical
        # to the previously-hardcoded version.
        extra_opts = " ".join(filter(None, [
            f"-E {config['tools']['dcc']['end_tol']}" if config["tools"]["dcc"]["end_tol"] != 5 else "",
            f"-m {config['tools']['dcc']['max_length']}" if config["tools"]["dcc"]["max_length"] != 1_000_000 else "",
            f"-n {config['tools']['dcc']['min_length']}" if config["tools"]["dcc"]["min_length"] != 30 else "",
        ])),
    threads: config["tools"]["dcc"]["threads"]
    log: "logs/dcc/{srr}.log"
    shell:
        """
        mkdir -p {params.outdir}
        ( cd {params.outdir} && DCC {input.junc} \
            -mt1 {input.junc1} \
            -mt2 {input.junc2} \
            -D -an {input.gtf} \
            -Pi {params.filter_f} {params.filter_m} -Nr {params.nr} \
            {params.extra_opts} \
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
        qc_bsj_threshold   = config["consensus"].get("qc_bsj_threshold", 5),
        adaptive_flag      = "--adaptive" if config["consensus"].get("adaptive", True) else "",
        adaptive_ratio     = config["consensus"].get("adaptive_ratio", 0.1),
        # Same key as the DCC -Nr first argument (rule dcc): the count given to CircCoordinates loci.
        dcc_count_threshold = config["tools"]["dcc"]["count_threshold"],
        # Only emit --min-confidence when it differs from consensus_filter.py's
        # own default (0.0), to keep the default-config command byte-identical.
        min_conf_arg       = (f"--min-confidence {config['consensus'].get('min_confidence', 0.0)}"
                               if config["consensus"].get("min_confidence", 0.0) != 0.0 else ""),
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
            --qc-bsj-threshold    {params.qc_bsj_threshold} \
            {params.adaptive_flag} \
            --adaptive-ratio      {params.adaptive_ratio} \
            --dcc-count-threshold {params.dcc_count_threshold} \
            {params.min_conf_arg} \
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
