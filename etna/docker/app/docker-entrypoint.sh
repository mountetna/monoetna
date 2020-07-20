#!/bin/bash
set -e

: "${RANDOM_MAX:=500}"
[ -n "$DEBUG" ] && echo "Running: $@"
export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"
cp npmrc ~/.npmrc || true

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  cd packages/etna-js
  npm install
fi

exec "$@"
