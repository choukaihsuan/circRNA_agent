"""
QC rules: FastQC on raw reads → fastp adapter trimming → MultiQC summary.
Tool paths resolved at DAG-build time via _find_tool() (defined in download.smk,
included before this file) so the rules work whether or not `conda activate`
was run before invoking Snakemake.
"""

_FASTQC  = _find_tool("fastqc")
_FASTP   = _find_tool("fastp")
_MULTIQC = _find_tool("multiqc")


rule fastqc_raw:
    input:
        r1 = RAW_DIR + "/{srr}_1.fastq.gz",
        r2 = RAW_DIR + "/{srr}_2.fastq.gz",
    output:
        html1 = RESULTS_DIR + "/qc/raw/{srr}_1_fastqc.html",
        zip1  = RESULTS_DIR + "/qc/raw/{srr}_1_fastqc.zip",
        html2 = RESULTS_DIR + "/qc/raw/{srr}_2_fastqc.html",
        zip2  = RESULTS_DIR + "/qc/raw/{srr}_2_fastqc.zip",
    params:
        outdir = RESULTS_DIR + "/qc/raw",
        fastqc = _FASTQC,
    threads: 2
    log: "logs/fastqc_raw/{srr}.log"
    shell:
        "{params.fastqc} {input.r1} {input.r2} -o {params.outdir} -t {threads} > {log} 2>&1"


rule fastp_trim:
    """Adapter trimming, quality filtering, and poly-G tail correction."""
    input:
        r1 = RAW_DIR + "/{srr}_1.fastq.gz",
        r2 = RAW_DIR + "/{srr}_2.fastq.gz",
    output:
        r1   = TRIMMED_DIR + "/{srr}_1.fastq.gz",
        r2   = TRIMMED_DIR + "/{srr}_2.fastq.gz",
        json = RESULTS_DIR + "/qc/fastp/{srr}.json",
        html = RESULTS_DIR + "/qc/fastp/{srr}.html",
    params:
        fastp = _FASTP,
    threads: config["threads"]
    log: "logs/fastp/{srr}.log"
    shell:
        """
        {params.fastp} \
            -i {input.r1}  -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -j {output.json} -h {output.html} \
            --thread {threads} \
            --detect_adapter_for_pe \
            --correction \
            --overrepresentation_analysis \
            --length_required 50 \
            --qualified_quality_phred 20 \
            > {log} 2>&1
        """


rule multiqc:
    """Aggregate FastQC and fastp reports into a single MultiQC report."""
    input:
        fastqc = expand(RESULTS_DIR + "/qc/raw/{srr}_1_fastqc.zip", srr=SAMPLES),
        fastp  = expand(RESULTS_DIR + "/qc/fastp/{srr}.json",        srr=SAMPLES),
    output:
        RESULTS_DIR + "/qc/multiqc_report.html",
    params:
        qc_dir  = RESULTS_DIR + "/qc",
        multiqc = _MULTIQC,
    log: "logs/multiqc.log"
    shell:
        "{params.multiqc} {params.qc_dir} -o {params.qc_dir} --filename multiqc_report.html -f > {log} 2>&1"
