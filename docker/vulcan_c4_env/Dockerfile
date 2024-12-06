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

# Install Snakemake
RUN pip install snakemake==8.10.4
RUN pip install snakemake-executor-plugin-slurm==0.4.2

# Download, configure, and install Slurm 21.08.4
RUN wget https://download.schedmd.com/slurm/slurm-21.08.8.tar.bz2 && \
    tar -xjf slurm-21.08.8.tar.bz2 && \
    cd slurm-21.08.8 && \
    ./configure --with-mysql_config=/usr/bin && \
    make && \
    make install && \
    cd .. && \
    rm -rf slurm-21.08.8.tar.bz2 slurm-21.08.8

# Add slurm dbd config file
COPY slurmdbd.conf /usr/local/etc/slurmdbd.conf
RUN chown root /usr/local/etc/slurmdbd.conf
RUN chmod 600 /usr/local/etc/slurmdbd.conf
