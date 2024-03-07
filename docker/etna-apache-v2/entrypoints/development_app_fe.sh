#!/usr/bin/env sh
set -e
set -x

if [ -z "$RELEASE_TEST" ]; then
  /entrypoints/build_app_fe.sh
fi

exec $@
