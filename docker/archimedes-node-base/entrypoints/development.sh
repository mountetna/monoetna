#!/usr/bin/env sh
set -e
set -x

if [ -z "$RELEASE_TEST" ]; then
  /entrypoints/build.sh
fi

if [ -e /app/development ]; then
  for hook in /app/development/*; do
    [ -x "$hook" ] && $hook
  done
fi

exec $@
