#!/usr/bin/env bash

set -e

# Fix SSH permissions after build.sh runs
echo "Copying SSH keys to root/.ssh"
cp -R /etc/ssh-keys-ro /root/.ssh 

echo "Setting permissions for root/.ssh"
chmod 700 /root/.ssh

[ -f /root/.ssh/config ] && chmod 600 /root/.ssh/config || true
[ -f /root/.ssh/known_hosts ] && chmod 644 /root/.ssh/known_hosts || true 

# Find all public key files and set permissions to 644 (rw-r--r--)
find /root/.ssh -name "*.pub" -type f -exec chmod 644 {} \;

# Find all private key files (no .pub extension) and set permissions to 600 (rw-------)
find /root/.ssh -type f ! -name "*.pub" ! -name "known_hosts*" ! -name "config" -exec chmod 600 {} \;

