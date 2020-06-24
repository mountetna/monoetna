#!/bin/bash
set -e

export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
  mkdir -p tmp/pids
  mkdir -p data/uploads
  mkdir -p data/data_blocks
  mkdir -p data/blueprints
  mkdir -p spec/data/stubs/blueprints
  mkdir -p spec/data/data_blocks

  rm -f tmp/pids/*.pid
  if [ -z "$SKIP_DB_WAIT" ]; then
    dockerize -wait tcp://metis_db:5432 -timeout 60s
    ./bin/metis migrate
  fi
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  npm --version
  node --version
  npm install --unsafe-perm
fi

exec "$@"
