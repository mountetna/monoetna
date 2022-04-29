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

mkdir /app/node_modules
chown -R 1000:1000 /app/node_modules

exec $@
