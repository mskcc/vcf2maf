FROM debian:bookworm-slim AS builder

# Install pacstrap and create a minimal base install in /install_root
RUN apt-get update && \
    apt-get install -y debootstrap ca-certificates curl && \
    debootstrap --variant=minbase --include=ca-certificates bookworm /install_root

# Download and install conda into the builder stage's /opt/conda
ENV MINICONDA_VERSION=py312_25.5.1-1
RUN curl -sL https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -bup /opt/conda && \
    rm -f /tmp/miniconda.sh

# Modify PATH to see conda, use faster dependency solver, and accept TOS
ENV PATH="/opt/conda/bin:${PATH}"
RUN conda config --set solver libmamba && \
    conda tos accept --channel defaults

# Use conda to install vcf2maf tools/dependencies into /usr/local
ENV VEP_VERSION=112.0 \
    HTSLIB_VERSION=1.20 \
    BCFTOOLS_VERSION=1.20 \
    SAMTOOLS_VERSION=1.20 \
    LIFTOVER_VERSION=447
RUN conda create -y -p /usr/local && \
    conda install -y -p /usr/local \
    -c conda-forge \
    -c bioconda \
    -c defaults \
    ensembl-vep==${VEP_VERSION} \
    htslib==${HTSLIB_VERSION} \
    bcftools==${BCFTOOLS_VERSION} \
    samtools==${SAMTOOLS_VERSION} \
    ucsc-liftover==${LIFTOVER_VERSION}

# Deploy the minimal OS and tools into a clean target layer
FROM scratch

LABEL maintainer="Cyriac Kandoth <ckandoth@gmail.com>"

COPY --from=builder /install_root /
COPY --from=builder /usr/local /usr/local
COPY data /opt/data
COPY *.pl /opt/
WORKDIR /opt
