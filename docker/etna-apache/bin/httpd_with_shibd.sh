#!/usr/bin/env bash

set -e
# Starts shibd in the background
/usr/sbin/shibd || true

# Start httpd
exec httpd -DFOREGROUND