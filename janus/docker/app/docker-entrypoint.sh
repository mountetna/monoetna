#!/bin/bash
set -e

: "${RANDOM_MAX:=500}"
[ -n "$DEBUG" ] && echo "Running: $@"
export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
  mkdir -p tmp/pids
  rm -f tmp/pids/*.pid
  if [ -z "$SKIP_DB_WAIT" ]; then
    dockerize -wait tcp://janus_db:5432 -timeout 60s
    ./bin/janus migrate
  fi
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  npm install --unsafe-perm
fi

exec "$@"
