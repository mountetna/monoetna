#!/usr/bin/env sh

set -e

export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"
# Default directories
mkdir -p /app/tmp/pids
mkdir -p /app/public/js
mkdir -p /app/public/css
mkdir -p /app/public/images
mkdir -p /app/log
mkdir -p /app/vendor/bundle
mkdir -p /app/data


if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
fi


if [ -n "$RUN_NPM_INSTALL" ]; then
  npm install --unsafe-perm
fi

if [ -e /app/build ]; then
  for hook in /app/build/*; do
    [ -e "$hook" ] && . "$hook"
  done
fi

cp /opt/npmrc /app/.npmrc
