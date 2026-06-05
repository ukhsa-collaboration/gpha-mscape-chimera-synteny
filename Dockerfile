FROM quay.io/condaforge/miniforge3:latest

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1

WORKDIR /app

RUN conda install -n base -c conda-forge -c bioconda -y python=3.12 git samtools s3cmd && conda clean --all --force-pkgs-dirs -y

COPY pyproject.toml README.md ./
COPY src ./src

RUN pip install .
