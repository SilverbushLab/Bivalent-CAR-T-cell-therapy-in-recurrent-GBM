# This file contains first Seurat and dependencies. Then installs infercnv and dependencies. The main issues to watch outfor is the rjags/jags/r-cran-rjags dependencies which are needed for infercnv.

# Next it installs python, jupyterlab, and necessary dependencies for scrublet.

# Use the Rocker RStudio base image. you have to use 4.3. 4.4 has broken signatures and will fail due to security for many packages.

FROM rocker/rstudio:4.3.3

# FIRST PART IS FOR SEURAT
#ALSO BIOCONDUCTOR AND TIDYVERSE

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE


# Install system dependencies for Seurat, FIt-SNE, and other packages

# Notes on inferCNV dependencies
# curl and below packages are infercnv dependencies. 
# r-cran-rjags IS KEY HERE.
# do not install jags from source and do not use libjags-dev!

RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    python3-venv \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    pkg-config \
    libgeos-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    curl \
    jags \
    r-cran-rjags \
    time \
    libgit2-dev \
    cmake \
    libssh2-1-dev \
    gfortran && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --- Python, venv, and Scrublet setup ---
# Create Python virtual environment and install dependencies. This is a way to do all installation into the virtual environment as this is required.
RUN python3 -m venv /venv && \
    /venv/bin/pip install --upgrade pip && \
    LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config && \
    /venv/bin/pip install llvmlite umap-learn numpy && \
    /venv/bin/pip install notebook ipykernel jupyterlab scrublet scipy scikit-learn matplotlib pandas anndata celltypist && \
    /venv/bin/pip install scanpy

ENV PATH="/venv/bin:${PATH}"

# --- Install R packages ---

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git && \
    g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp -o /usr/local/bin/fast_tsne -pthread -lfftw3 -lm && \
    rm -rf FIt-SNE

# Install R and Bioconductor dependencies. Warning this takes a long time.
RUN R --no-save -e "install.packages('BiocManager')"

RUN R --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi'))" 

RUN R --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools', 'spatstat.explore', 'spatstat.geom', 'hdf5r', 'rgeos', 'remotes', 'tidyverse', 'readr'))"

RUN R --no-save -e "install.packages('Matrix')"

RUN R --no-save -e "install.packages('Seurat')"

#dependencies for scCustomize
RUN apt-get update && apt-get install -y \
    libcairo2-dev \
    libjpeg-dev \
    libtiff5-dev \
    libfreetype6-dev

RUN R -e "install.packages('scCustomize')"

RUN R -e 'remotes::install_github("mojaveazure/seurat-disk")'

# new packages
RUN R -e "install.packages(c('devtools', 'patchwork', 'gert', 'harmony', 'gghalves', 'ggpubr'))"

RUN R -e "devtools::install_github('immunogenomics/presto')"

# -----InferCNV install -----


# Install additional R dependencies required for infercnv
RUN R --no-echo --no-restore --no-save -e "install.packages(c('phyclust', 'fastcluster', 'parallelDist', 'HiddenMarkov', 'foreach', 'doParallel', 'coda', 'argparse', 'coin'), repos = 'http://cran.us.r-project.org')" && \
    R --no-echo --no-restore --no-save -e "BiocManager::install('edgeR')"
RUN R --no-echo --no-restore --no-save -e "install.packages(c('rjags'), repos = 'http://cran.us.r-project.org')" # last step is to install the rjags pagacke in R.


# Clone and install infercnv
RUN git clone https://github.com/broadinstitute/infercnv && \
    cd infercnv && \
    git checkout master && \
    git checkout 1b46b48303bac4a882bcb758e78fcf7f832fdefb && \
    R CMD INSTALL . && \
    rm -rf infercnv
 

ENV PATH=${PATH}:/infercnv/scripts

# -----scRepertoire install -----

# Need reticulate for R
RUN R -e "install.packages('reticulate', repos = 'https://cloud.r-project.org')"

#-- environment setup --

# Download and install Miniforge (aarch64 version)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:${PATH}"

# Create r-reticulate conda environment with Python 3.9
RUN conda create -y -n r-reticulate python=3.9

# (optional) Set ownership to RStudio user if needed
RUN chown -R rstudio /opt/conda

# complex setup of dependencies to enable latest version of scRepertoire. In dockerfile will need to setup r-reticulat env before doing this. Use miniforge to do this.

RUN conda run -n r-reticulate pip install keras tensorflow "jax[cpu]" \
    umap-learn numpy notebook ipykernel jupyterlab scrublet scipy scikit-learn \
    matplotlib pandas anndata celltypist scanpy

# setup and install

# Step 1: Load reticulate and force the correct Python env
RUN R -e "Sys.setenv(PATH = paste('/opt/conda/bin', Sys.getenv('PATH'), sep = ':')); \
          library(reticulate); \
          use_condaenv('r-reticulate', required = TRUE)"


# Step 2: install dependencies for immApex
#RUN R -e "py_install(c('keras', 'tensorflow', 'jax[cpu]''), pip = TRUE)"  # needed dependencies for immApex

# Step 3: install immApex
RUN R -e "devtools::install_github('BorchLab/immApex')"

#Step 4: install scRepertoire
RUN R -e "devtools::install_github('BorchLab/scRepertoire', ref = 'v2.3.4')"

