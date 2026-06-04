# Use a base image with Python installed, as Snakemake is a Python package
FROM python:3.11-slim

# Install SSH server, build dependencies, core utilities, Slurm dependencies, and MariaDB
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        openssh-server \
        build-essential \
        libffi-dev \
        libssl-dev \
        coreutils \
        git \
        munge \
        wget \
        libmunge-dev \
        libmunge2 \
        libreadline-dev \
        zlib1g-dev \
        mariadb-server \
        mariadb-client \
        libmariadb-dev && \
    echo 'root:root' | chpasswd

# Setup SSH server
RUN mkdir /var/run/sshd
RUN sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# Expose the default SSH port
EXPOSE 22

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /$HOME/miniconda && \
    rm miniconda.sh

# Create conda environmnent
RUN . /$HOME/miniconda/etc/profile.d/conda.sh && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    conda create -n conda-vulcan -c conda-forge -c bioconda snakemake=8.10.4 snakemake-executor-plugin-slurm=0.4.2 -y

# TODO inject conda func into bashrc

# Download, configure, and install Slurm 21.08.4
RUN wget https://download.schedmd.com/slurm/slurm-21.08.8.tar.bz2 && \
    tar -xjf slurm-21.08.8.tar.bz2 && \
    cd slurm-21.08.8 && \
    ./configure --with-mysql_config=/usr/bin && \
    make && \
    make install && \
    cd .. && \
    rm -rf slurm-21.08.8.tar.bz2 slurm-21.08.8
