#!/bin/bash
set -e

export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
  mkdir -p tmp/pids

  rm -f tmp/pids/*.pid
fi

exec "$@"
