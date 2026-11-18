# syntax=docker/dockerfile:1.6

FROM continuumio/miniconda3:latest

LABEL maintainer="Brown Beckley <brownbeckley94@gmail.com>"
LABEL description="Kleboscope - Comprehensive Klebsiella pneumoniae Genomic Typing Platform"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        jq \
        git \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/kleboscope

COPY environment.yml /opt/kleboscope/environment.yml

RUN --mount=type=cache,target=/opt/conda/pkgs \
    conda env create -f /opt/kleboscope/environment.yml

RUN conda install -y -n kleboscope -c conda-forge -c bioconda fastani

ENV CONDA_PREFIX=/opt/conda/envs/kleboscope
ENV CONDA_DEFAULT_ENV=kleboscope
ENV PATH=/opt/conda/envs/kleboscope/bin:${PATH}

COPY . /opt/kleboscope/

RUN pip install --no-deps --no-build-isolation -e /opt/kleboscope

RUN chmod -R a+rX /opt/kleboscope/kleboscope/modules/ && \
    find /opt/kleboscope/kleboscope/modules/ -type f \
        \( -name "*.py" -o -name "*.pl" -o -name "*.sh" \) \
        -exec chmod a+x {} \;

RUN abricate --setupdb

RUN test -d /opt/kleboscope/kleboscope/modules/mlst_module/db/pubmlst/klebsiella || \
    { echo "[ERROR] MLST scheme missing from build context (kleboscope/modules/mlst_module/db/)"; exit 1; } && \
    ls /opt/kleboscope/kleboscope/modules/mlst_module/db/pubmlst/klebsiella/*.tfa >/dev/null 2>&1 || \
    { echo "[ERROR] No allele files (*.tfa) in packaged MLST scheme"; exit 1; }

RUN kleboscope --pull-mlst-db && \
    test -f /opt/kleboscope/kleboscope/modules/mlst_module/db/blast/mlst.fa || \
    { echo "[ERROR] MLST blast index missing after index build — check makeblastdb in environment.yml"; exit 1; }

RUN test -f /opt/kleboscope/kleboscope/modules/mge_module/mobileOG-db/beatrix-1-6_v1_all/mobileOG-db-beatrix-1.6.dmnd || \
    { echo "[ERROR] mobileOG-db DIAMOND index missing from build context (kleboscope/modules/mge_module/mobileOG-db/)"; exit 1; }

RUN kleboscope --pull-mge-db

RUN groupadd -g 1000 kleboscope 2>/dev/null || true && \
    useradd -u 1000 -g 1000 -m -s /bin/bash kleboscope 2>/dev/null || true && \
    chown -R 1000:1000 /opt/kleboscope

ENV HOME=/tmp
ENV MAMBA_ROOT_PREFIX=/opt/conda

RUN printf '%s\n' \
    '#!/bin/bash' \
    'set -e' \
    'export CONDA_PREFIX=/opt/conda/envs/kleboscope' \
    'export CONDA_DEFAULT_ENV=kleboscope' \
    'export PATH=${CONDA_PREFIX}/bin:${PATH}' \
    'export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}' \
    'if [ $# -gt 0 ] && [ "${1#-}" != "$1" ]; then' \
    '    set -- kleboscope "$@"' \
    'fi' \
    'exec "$@"' \
    > /usr/local/bin/entrypoint.sh \
    && chmod +x /usr/local/bin/entrypoint.sh

USER 1000:1000

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["kleboscope", "-h"]
