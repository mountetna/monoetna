# Use a base image with Python installed, as Snakemake is a Python package
FROM python:3.11-slim

# Install SSH server, build dependencies, and core utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        openssh-server \
        build-essential \
        libffi-dev \
        libssl-dev \
        coreutils \
        git && \
    echo 'root:root' | chpasswd

# Setup SSH server
RUN mkdir /var/run/sshd
RUN sed -i 's/#PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# Expose the default SSH port
EXPOSE 22

# Install Snakemake
RUN pip install snakemake==8.4.0

# Start SSH server as the default command
CMD ["/usr/sbin/sshd", "-D"]
