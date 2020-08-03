#!/usr/bin/env sh
set -e

/app/entrypoints/build.sh

export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

rm -f tmp/pids/*.pid

if [ -z "$SKIP_RUBY_SETUP" ]; then

  if [ -z "$SKIP_DB_WAIT" ]; then
    dockerize -wait tcp://${APP_NAME}_db:5432 -timeout 60s
    ./bin/${APP_NAME} migrate
  fi
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  npm link ../etna/packages/etna-js
fi

if [ -e /app/development ]; then
  for hook in /app/development/*; do
    [ -e "$hook" ] && $hook
  done
fi

exec "$@"
