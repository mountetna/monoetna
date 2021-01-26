#!/usr/bin/env bash

set -e
# Starts shibd in the background
/usr/sbin/shibd

# Start httpd
exec httpd -DFOREGROUND