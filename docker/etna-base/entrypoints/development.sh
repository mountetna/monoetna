#!/usr/bin/env sh
set -e
set -x

if [ -z "$SKIP_BUILD" ]; then
  if [ -z "$RELEASE_TEST" ]; then
    /entrypoints/build.sh
  fi

  export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"

  rm -f tmp/pids/*.pid

  if [ -e /app/development ]; then
    for hook in /app/development/*; do
      [ -x "$hook" ] && $hook
    done
  fi

  if [ -z "$SKIP_RUBY_SETUP" ]; then
    if [ -z "$SKIP_DB_WAIT" ]; then
      dockerize -wait tcp://${APP_NAME}_db:5432 -timeout 60s
      ./bin/${APP_NAME} migrate
    fi
  fi

  mkdir -p /app/data/uploads
  mkdir -p /app/data/data_blocks
fi

if [ -n "$WAIT_FOR_APP" ]; then
  dockerize -wait tcp://${APP_NAME}_app:3000 -timeout 60s
fi

exec $@
