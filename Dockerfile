FROM clearlinux:latest AS builder

ENV CLEAR_VERSION=32820 \
    MINICONDA_VERSION=py38_4.8.2 \
    VEP_VERSION=99.2

# Prepare a root directory with minimal OS and some dependencies
RUN swupd os-install --no-progress --no-boot-update --no-scripts \
    --version ${CLEAR_VERSION} \
    --path /install_root \
    --statedir /swupd-state \
    --bundles os-core-update,which,samtools

# Use conda to install VEP and remaining dependencies into /usr/local
RUN swupd bundle-add --no-progress curl && \
    curl -sL https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    sh /tmp/miniconda.sh -bfp /usr && \
    conda create -qy -p /usr/local -c conda-forge -c bioconda -c defaults ensembl-vep==${VEP_VERSION}

# Deploy the minimal OS and tools into a clean target layer
FROM scratch

LABEL maintainer="Cyriac Kandoth <ckandoth@gmail.com>"

COPY --from=builder /install_root /
COPY --from=builder /usr/local /usr/local
COPY data /opt/data
COPY *.pl /opt/
WORKDIR /opt
