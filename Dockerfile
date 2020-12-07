FROM clearlinux:latest AS builder

# Install a minimal versioned OS into /install_root, and bundled tools if any
ENV CLEAR_VERSION=33980
RUN swupd os-install --no-progress --no-boot-update --no-scripts \
    --version ${CLEAR_VERSION} \
    --path /install_root \
    --statedir /swupd-state \
    --bundles os-core-update,which

# Download and install conda into /usr/bin
ENV MINICONDA_VERSION=py37_4.9.2
RUN swupd bundle-add --no-progress curl && \
    curl -sL https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    sh /tmp/miniconda.sh -bfp /usr

# Use conda to install remaining tools/dependencies into /usr/local
ENV VEP_VERSION=102.0 \
    HTSLIB_VERSION=1.10.2 \
    BCFTOOLS_VERSION=1.10.2 \
    SAMTOOLS_VERSION=1.10 \
    LIFTOVER_VERSION=377
RUN conda create -qy -p /usr/local \
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
