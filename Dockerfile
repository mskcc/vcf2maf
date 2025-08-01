FROM cachyos/cachyos:latest

# Download and install conda into /usr/local
ENV MINICONDA_VERSION=py312_25.5.1-1
RUN curl -sL https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -bup /usr/local && \
    rm -f /tmp/miniconda.sh && \
    conda config --set solver libmamba && \
    conda tos accept --channel defaults

# Use conda to install remaining tools/dependencies into /usr/local
ENV VEP_VERSION=112.0 \
    HTSLIB_VERSION=1.20 \
    BCFTOOLS_VERSION=1.20 \
    SAMTOOLS_VERSION=1.20 \
    LIFTOVER_VERSION=447
RUN conda install -y -p /usr/local \
    -c conda-forge \
    -c bioconda \
    -c defaults \
    ensembl-vep==${VEP_VERSION} \
    htslib==${HTSLIB_VERSION} \
    bcftools==${BCFTOOLS_VERSION} \
    samtools==${SAMTOOLS_VERSION} \
    ucsc-liftover==${LIFTOVER_VERSION}

LABEL maintainer="Cyriac Kandoth <ckandoth@gmail.com>"

COPY data /opt/data
COPY *.pl /opt/
WORKDIR /opt
