#!/bin/bash

# Initialize the Git repository if it doesn't exist
if [ ! -d "/app/available-workflows/test-repo/.git" ]; then
    git config --global user.email "you@example.com"
    git config --global user.name "Your Name"
    git config --global init.defaultBranch main
    git init /app/available-workflows/test-repo
    git -C /app/available-workflows/test-repo add .
    git -C /app/available-workflows/test-repo commit -m "Initial commit"
fi

# Re-run sshd
exec /usr/sbin/sshd -D