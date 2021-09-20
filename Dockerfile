FROM docker.io/ubuntu:20.04

LABEL maintaner="Anže Lovše (@AnzeLovse)"
LABEL release-date="2021-06-22"

RUN echo "debconf debconf/frontend select Noninteractive" | debconf-set-selections \
    && apt-get update  \
    && apt-get -qq -y install \
    apt-utils \
    curl \
    bzip2 \
    software-properties-common \
    libcurl4-openssl-dev \
    build-essential \
    libssl-dev \
    libxml2-dev

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
    && apt-get -qq -y install r-base r-base-dev

RUN R -e "install.packages(c('BiocManager'), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c('DESeq2', 'apeglm'))"
RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh

COPY environment.yaml /var/cache/build/
RUN conda update -y conda \
    && conda install -y -c conda-forge mamba \
    && mamba create -q -y -c conda-forge -c bioconda -n snakemake snakemake python=3.8 \
    && conda env update -n snakemake --file /var/cache/build/environment.yaml \
    && conda run -n snakemake python -m pip install rnanorm

RUN apt-get clean \
    && conda clean --all -y \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && rm -rf /tmp/*
