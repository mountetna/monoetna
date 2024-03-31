# Use a base image with Python installed, as Snakemake is a Python package
FROM python:3.11-slim

# Install SSH server, build dependencies, core utilities, and Slurm dependencies
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
        zlib1g-dev && \
    echo 'root:root' | chpasswd

# Setup SSH server
RUN mkdir /var/run/sshd
RUN sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# Expose the default SSH port
EXPOSE 22

# Install Snakemake
RUN pip install snakemake==8.4.0

# Download, configure, and install Slurm 21.08.4
RUN wget https://download.schedmd.com/slurm/slurm-21.08.8.tar.bz2 && \
    tar -xjf slurm-21.08.8.tar.bz2 && \
    cd slurm-21.08.8 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf slurm-21.08.8.tar.bz2 slurm-21.08.8

# Start SSH server as the default command
CMD ["/usr/sbin/sshd", "-D"]
