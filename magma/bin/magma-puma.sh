#!/usr/bin/env bash

# Magma restarts Puma after some model updates. In production, that restart can
# leave Puma unavailable while the container is still otherwise healthy.
#
# This wrapper runs as PID 1 and starts Puma as a child process. If Puma exits,
# the wrapper starts a fresh Puma process without requiring a full container
# restart. If Kubernetes stops the container, the wrapper forwards that shutdown
# signal to Puma and exits.

set -u

PUMA_CONFIG="/entrypoints/puma.rb"
RESTART_DELAY_SECONDS=2

puma_pid=""

shutdown() {
  if [ -n "$puma_pid" ]; then
    echo "Stopping Puma pid ${puma_pid}"
    kill -TERM "$puma_pid" 2>/dev/null
    wait "$puma_pid" 2>/dev/null
  fi

  exit 0
}

# Kubernetes stops the container by sending SIGTERM to PID 1; forward it to Puma.
trap shutdown TERM INT

while true; do
  # Run Puma as a child so this wrapper can restart it without restarting the container.
  echo "Starting Puma"
  puma -C "$PUMA_CONFIG" &
  puma_pid=$!

  wait "$puma_pid"
  puma_status=$?
  puma_pid=""

  # The loop repeats after this pause, starting a fresh Puma process.
  echo "Puma exited with status ${puma_status}; restarting in ${RESTART_DELAY_SECONDS}s"
  sleep "$RESTART_DELAY_SECONDS"
done
