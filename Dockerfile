FROM clearlinux:latest AS builder

# Install a minimal versioned OS into /install_root, and bundled tools if any
ENV CLEAR_VERSION=41780
RUN swupd os-install --no-progress --no-boot-update --no-scripts \
    --version ${CLEAR_VERSION} \
    --path /install_root \
    --statedir /swupd-state \
    --bundles os-core-update,which

# Download and install conda into /usr/bin
ENV MINICONDA_VERSION=py312_24.4.0-0
RUN curl -sL https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -bup /usr && \
    rm -f /tmp/miniconda.sh && \
    conda config --set solver libmamba

# Use mamba to install remaining tools/dependencies into /usr/local
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
