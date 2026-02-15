FROM continuumio/miniconda3

WORKDIR /app

RUN apt-get update \
    # install linux dependencies
    && apt-get install -y --no-install-recommends wget bash build-essential pip unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    # install PLINK
    && wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250707.zip \
    && unzip plink2_linux_x86_64_20250707.zip \
    && rm -rf plink2_linux_x86_64_20250707.zip \
    # move plink2 executable to $PATH
    && mv plink2 /usr/bin \
    # install fraposa- automatically on path
    && pip install fraposa \
    # install pgs catalog ancestry package- automatically on path
    && pip install pgscatalog-ancestry-adjust \
    # install python and python packages
    && conda install -y -n base -c conda-forge -c bioconda scipy pandas python=3.11 \
    && conda clean --all --yes \

USER root