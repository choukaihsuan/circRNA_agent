"""
QC rules: FastQC on raw reads → fastp adapter trimming → MultiQC summary.
"""


rule fastqc_raw:
    input:
        r1 = RAW_DIR + "/{srr}_1.fastq.gz",
        r2 = RAW_DIR + "/{srr}_2.fastq.gz",
    output:
        html1 = "results/qc/raw/{srr}_1_fastqc.html",
        zip1  = "results/qc/raw/{srr}_1_fastqc.zip",
        html2 = "results/qc/raw/{srr}_2_fastqc.html",
        zip2  = "results/qc/raw/{srr}_2_fastqc.zip",
    threads: 2
    log: "logs/fastqc_raw/{srr}.log"
    shell:
        "fastqc {input.r1} {input.r2} -o results/qc/raw/ -t {threads} > {log} 2>&1"


rule fastp_trim:
    """Adapter trimming, quality filtering, and poly-G tail correction."""
    input:
        r1 = RAW_DIR + "/{srr}_1.fastq.gz",
        r2 = RAW_DIR + "/{srr}_2.fastq.gz",
    output:
        r1   = TRIMMED_DIR + "/{srr}_1.fastq.gz",
        r2   = TRIMMED_DIR + "/{srr}_2.fastq.gz",
        json = "results/qc/fastp/{srr}.json",
        html = "results/qc/fastp/{srr}.html",
    threads: config["threads"]
    log: "logs/fastp/{srr}.log"
    shell:
        """
        fastp \
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
        fastqc = expand("results/qc/raw/{srr}_1_fastqc.zip", srr=SAMPLES),
        fastp  = expand("results/qc/fastp/{srr}.json",        srr=SAMPLES),
    output:
        "results/qc/multiqc_report.html",
    log: "logs/multiqc.log"
    shell:
        "multiqc results/qc/ -o results/qc/ --filename multiqc_report.html -f > {log} 2>&1"
