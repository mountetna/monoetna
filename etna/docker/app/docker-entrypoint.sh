#!/bin/bash
set -e

: "${RANDOM_MAX:=500}"
[ -n "$DEBUG" ] && echo "Running: $@"
export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  npm install --unsafe-perm
fi

exec "$@"
