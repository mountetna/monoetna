#!/bin/bash

# Initialize the Git repository if it doesn't exist
if [ ! -d "/test-utils/available-workflows/snakemake-repo/.git" ]; then
    git config --global user.email "you@example.com"
    git config --global user.name "Your Name"
    git config --global init.defaultBranch main
    git init /test-utils/available-workflows/snakemake-repo
    git -C /test-utils/available-workflows/snakemake-repo add .
    git -C /test-utils/available-workflows/snakemake-repo commit -m "Initial commit"
    git -C /test-utils/available-workflows/snakemake-repo tag v1
fi