#!/usr/bin/env sh
set -e

#if [ -n "$VERBOSE" ]; then
  set -x
#fi

if [ -z "$SKIP_BUILD" ]; then
  if [ -z "$RELEASE_TEST" ]; then
    /entrypoints/build.sh
  fi

  if [ -e /app/development ]; then
    for hook in /app/development/*; do
      [ -x "$hook" ] && $hook
    done
  fi
fi

exec poetry run $@
