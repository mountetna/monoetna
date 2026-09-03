#!/usr/bin/env bash

# Watch the puma.hold file for a single Puma process. If the hold file expires,
# assume Puma did not finish restarting and kill it so magma-puma.sh can start a
# fresh process.

set -u

if [ "$#" -lt 1 ]; then
  echo "Usage: puma-hold-monitor.sh <puma-pid>"
  exit 1
fi

PUMA_PID="$1"
HOLD_FILE=$(grep -E "^\s+:hold_file:" config.yml | awk '{ print $2 }')
MONITOR_INTERVAL_SECONDS=1

[ -z "$HOLD_FILE" ] && exit 0

while true; do
  if ! kill -0 "$PUMA_PID" 2>/dev/null; then
    exit 0
  fi

  if [ -f "$HOLD_FILE" ]; then
    hold_time=$(cat "$HOLD_FILE")
    curr_time=$(date -u --iso-8601=sec)

    if [[ ! $curr_time < $hold_time ]]; then
      echo "puma.hold expired; killing Puma pid ${PUMA_PID}"
      kill -9 "$PUMA_PID" 2>/dev/null
      rm -f "$HOLD_FILE"
      exit 0
    fi
  fi

  sleep "$MONITOR_INTERVAL_SECONDS"
done
