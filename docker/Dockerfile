# Lightweight toolkit image for streaming NA12878 (or any CRAM) and calling ESR1/ESR2 without local alignment
#
# Usage (from repo root):
#   docker build -t ad-finemap -f docker/Dockerfile .
#   docker run --rm -it -v $PWD:/work -w /work ad-finemap \
#     bash scripts/stream_call.sh chr6:146809107-157129619 NA12878.ESR1
#
# Requirements on host: Docker/Podman and network access to the CRAM URL.

FROM mambaorg/micromamba:1.5.8

USER root
SHELL ["/bin/bash", "-c"]

# Minimal OS deps
RUN apt-get update && apt-get install -y --no-install-recommends curl unzip ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Install core tools from bioconda/conda-forge.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    samtools bcftools htslib plink \
    r-base r-data.table r-ggplot2 bioconductor-geoquery bioconductor-limma && \
    micromamba clean -a -y

WORKDIR /work
ENTRYPOINT ["/bin/bash"]
