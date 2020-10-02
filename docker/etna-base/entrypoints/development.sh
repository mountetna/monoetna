#!/usr/bin/env sh
set -e

#if [ -n "$VERBOSE" ]; then
  set -x
#fi

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
      #Sometimes during startup the cpu can be slammed for all the servies, so we need to wait
      # a decent amount of time for all pieces to come up.
      dockerize -wait tcp://${APP_NAME}_db:5432 -timeout 300s
      ./bin/${APP_NAME} migrate
    fi
  fi

  mkdir -p /app/data/uploads
  mkdir -p /app/data/data_blocks
fi

if [ -n "$WAIT_FOR_APP" ]; then
  #Sometimes during startup the cpu can be slammed for all the servies, so we need to wait
  # a decent amount of time for all pieces to come up.
  dockerize -wait tcp://${APP_NAME}_app:3000 -timeout 300s
fi

exec $@
