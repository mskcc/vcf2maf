FROM clearlinux:latest AS builder

# Install a minimal versioned OS into /install_root, and bundled tools if any
ENV CLEAR_VERSION=32820
RUN swupd os-install --no-progress --no-boot-update --no-scripts \
    --version ${CLEAR_VERSION} \
    --path /install_root \
    --statedir /swupd-state \
    --bundles os-core-update,which

# Download and install conda into /usr/bin
ENV MINICONDA_VERSION=py38_4.8.2
RUN swupd bundle-add --no-progress curl && \
    curl -sL https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    sh /tmp/miniconda.sh -bfp /usr

# Use conda to install remaining tools/dependencies into /usr/local
ENV VEP_VERSION=99.2 \
    HTSLIB_VERSION=1.9 \
    BCFTOOLS_VERSION=1.9 \
    SAMTOOLS_VERSION=1.9
RUN conda create -qy -p /usr/local \
        -c conda-forge \
        -c bioconda \
        -c defaults \
        ensembl-vep==${VEP_VERSION} \
        htslib==${HTSLIB_VERSION} \
        bcftools==${BCFTOOLS_VERSION} \
        samtools==${SAMTOOLS_VERSION}

# Deploy the minimal OS and tools into a clean target layer
FROM scratch

LABEL maintainer="Cyriac Kandoth <ckandoth@gmail.com>"

COPY --from=builder /install_root /
COPY --from=builder /usr/local /usr/local
COPY data /opt/data
COPY *.pl /opt/
WORKDIR /opt
