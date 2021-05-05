FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="16f31f1d9470da11b5ae066f139c71cda6aa5d94cdd1266e0de8c211c1040247"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: ../workflow/envs/featurecounts.yaml
#   prefix: /conda-envs/f970be3bdf66291746913c7fe60a78a4
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - subread ==2.0.1
RUN mkdir -p /conda-envs/f970be3bdf66291746913c7fe60a78a4
COPY ../workflow/envs/featurecounts.yaml /conda-envs/f970be3bdf66291746913c7fe60a78a4/environment.yaml

# Conda environment:
#   source: ../workflow/envs/samtools.yaml
#   prefix: /conda-envs/9608721699f97513ba7f47bd4e3db24b
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - samtools ==1.10
RUN mkdir -p /conda-envs/9608721699f97513ba7f47bd4e3db24b
COPY ../workflow/envs/samtools.yaml /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.72.0-5-gd0d054c/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/200f469099e32bb3d83d9ca89b039191
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bowtie2 ==2.4.1  # Keep consistent with version specified in bowtie2/align
#     - samtools ==1.10
RUN mkdir -p /conda-envs/200f469099e32bb3d83d9ca89b039191
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.72.0-5-gd0d054c/bio/bowtie2/build/environment.yaml /conda-envs/200f469099e32bb3d83d9ca89b039191/environment.yaml

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/06aa832cdf405d93f3eacda03259215e
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bowtie2 ==2.4.1  # Keep consistent with version specified in bowtie2/build
#     - samtools ==1.10
RUN mkdir -p /conda-envs/06aa832cdf405d93f3eacda03259215e
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/bowtie2/align/environment.yaml /conda-envs/06aa832cdf405d93f3eacda03259215e/environment.yaml

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/cutadapt/pe/environment.yaml
#   prefix: /conda-envs/473f10b47f8c5cd1ed7bd3701ad805a5
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - cutadapt ==2.10
RUN mkdir -p /conda-envs/473f10b47f8c5cd1ed7bd3701ad805a5
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/cutadapt/pe/environment.yaml /conda-envs/473f10b47f8c5cd1ed7bd3701ad805a5/environment.yaml

# Conda environment:
#   source: https:/github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/08d4368302a4bdf7eda6b536495efe7d
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - fastqc ==0.11.9
RUN mkdir -p /conda-envs/08d4368302a4bdf7eda6b536495efe7d
ADD https://github.com/snakemake/snakemake-wrappers/raw/0.72.0/bio/fastqc/environment.yaml /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/f970be3bdf66291746913c7fe60a78a4 --file /conda-envs/f970be3bdf66291746913c7fe60a78a4/environment.yaml && \
    mamba env create --prefix /conda-envs/9608721699f97513ba7f47bd4e3db24b --file /conda-envs/9608721699f97513ba7f47bd4e3db24b/environment.yaml && \
    mamba env create --prefix /conda-envs/200f469099e32bb3d83d9ca89b039191 --file /conda-envs/200f469099e32bb3d83d9ca89b039191/environment.yaml && \
    mamba env create --prefix /conda-envs/06aa832cdf405d93f3eacda03259215e --file /conda-envs/06aa832cdf405d93f3eacda03259215e/environment.yaml && \
    mamba env create --prefix /conda-envs/473f10b47f8c5cd1ed7bd3701ad805a5 --file /conda-envs/473f10b47f8c5cd1ed7bd3701ad805a5/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba clean --all -y
