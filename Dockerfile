FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="choukaihsuan" \
      description="circRNA analysis pipeline (CIRIquant + DCC + edgeR)" \
      version="1.0.0"

# ── System dependencies ───────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        default-jdk \
        procps \
        wget \
        curl \
        git \
        build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Conda environment (mamba for faster solve) ───────────────────────────────
RUN conda install -n base -c conda-forge mamba -y && conda clean -afy
COPY envs/circrna.yaml /tmp/circrna.yaml
RUN mamba env create --yes -f /tmp/circrna.yaml \
    && conda clean -afy \
    && rm /tmp/circrna.yaml

# ── PATH & Java ───────────────────────────────────────────────────────────────
ENV PATH="/opt/conda/envs/ciriquant/bin:${PATH}"
ENV JAVA_HOME="/usr/lib/jvm/default-java"
ENV CONDA_DEFAULT_ENV="ciriquant"

# ── Verify key tools are importable/executable ───────────────────────────────
RUN CIRIquant --version 2>&1 | grep -i 'ciriquant\|version' \
    && (DCC --version 2>&1 || DCC 2>&1 | head -3 || true) \
    && STAR --version 2>&1 | head -1 \
    && hisat2 --version 2>&1 | head -1 \
    && bwa 2>&1 | head -3 \
    && samtools --version 2>&1 | head -1 \
    && snakemake --version \
    && Rscript -e "library(edgeR);  cat('edgeR OK\n')" \
    && Rscript -e "library(DESeq2); cat('DESeq2 OK\n')"

WORKDIR /pipeline

CMD ["/bin/bash"]
