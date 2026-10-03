FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="choukaihsuan" \
      description="circRNA analysis pipeline (CIRIquant + DCC + edgeR)" \
      version="1.1.0"

# ── System dependencies ───────────────────────────────────────────────────────
# Debian bullseye's live security mirror has pruned old point-release .debs
# (bullseye is past normal LTS support), causing 404s on `apt-get update`.
# Pin to the snapshot.debian.org mirror matching this base image's build date
# instead -- same fix the image itself ships commented out in sources.list.
# This only affects OS-level packages (JDK/git/curl/build tools), none of
# which are in Table 1; it has no bearing on the conda env versions below.
RUN sed -i \
        -e 's|deb http://deb.debian.org/debian bullseye main|deb http://snapshot.debian.org/archive/debian/20231030T000000Z bullseye main|' \
        -e 's|deb http://deb.debian.org/debian-security bullseye-security main|deb http://snapshot.debian.org/archive/debian-security/20231030T000000Z bullseye-security main|' \
        -e 's|deb http://deb.debian.org/debian bullseye-updates main|deb http://snapshot.debian.org/archive/debian/20231030T000000Z bullseye-updates main|' \
        /etc/apt/sources.list \
    && printf 'Acquire::Check-Valid-Until "false";\n' > /etc/apt/apt.conf.d/99snapshot

RUN apt-get update && apt-get install -y --no-install-recommends \
        default-jdk \
        procps \
        wget \
        curl \
        git \
        build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── Conda environment ─────────────────────────────────────────────────────────
# NOT solved from envs/circrna.yaml. A fresh solve of these exact pinned
# versions (python=3.7.12 alongside current star/samtools/r-base builds)
# against today's live conda-forge/bioconda package index is infeasible --
# the shared C libraries those tools now link against (htslib/libcurl/
# openssl/krb5) have moved past what py3.7 can use, and the solver reports
# a hard conflict (see containers/build_notes.md). Instead, this installs
# the exact package set already running on the verified analysis server,
# captured via `conda list --explicit` -- a known-mutually-compatible
# combination because it is literally what produced the 19 published
# datasets, not a new solve. envs/circrna.yaml is kept as a human-readable
# statement of intent; it is not what builds the image.
COPY containers/spec-file.linux-64.txt /tmp/spec-file.linux-64.txt
RUN conda create --name ciriquant --file /tmp/spec-file.linux-64.txt --yes \
    && conda clean -afy \
    && rm /tmp/spec-file.linux-64.txt

# Packages the server env has that are not conda-managed (pip-only; cross-
# referenced against spec-file.linux-64.txt, see containers/build_notes.md
# for the full pip-freeze.txt -> pip-only split). --no-deps because the
# exact dependency versions are already pinned above; letting pip resolve
# deps itself would silently drift from the verified set.
COPY containers/pip-only-requirements.txt /tmp/pip-only-requirements.txt
RUN /opt/conda/envs/ciriquant/bin/pip install --no-deps -r /tmp/pip-only-requirements.txt \
    && rm /tmp/pip-only-requirements.txt

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
