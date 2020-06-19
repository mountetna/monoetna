#!/bin/bash
set -e

: "${RANDOM_MAX:=500}"
[ -n "$DEBUG" ] && echo "Running: $@"
export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
  mkdir -p tmp/pids
  mkdir -p spec/data/stubs/blueprints
  mkdir -p spec/data/data_blocks

  rm -f tmp/pids/*.pid
  if [ -z "$SKIP_DB_WAIT" ]; then
    dockerize -wait tcp://metis_db:5432 -timeout 60s
    ./bin/metis migrate
  fi
else
  while ! bundle check >/dev/null 2>&1; do
    echo "Awaiting for make bundle on host..."
    sleep 5
  done
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  npm install --unsafe-perm
fi

exec "$@"
