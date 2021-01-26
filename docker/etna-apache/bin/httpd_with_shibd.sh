#!/usr/bin/env bash

cleanup() {
  kill -- -$(jobs -p)
}

set -e
set -m

# Starts shibd in the background
/usr/sbin/shibd -f -F &>/tmp/shibd.log &
SHIBD_PID=$!
# Start httpd
httpd -DFOREGROUND &>/tmp/httpd.log &
HTTPD_PID=$!
trap cleanup EXIT
tail -f /tmp/*.log /var/log/shibboleth/shibd.log

