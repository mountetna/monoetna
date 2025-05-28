#!/usr/bin/env bash
set -e
set -x

app_name_capitalized=$(echo ${APP_NAME} | tr [a-z] [A-Z])

# In which case, we are using development.sh as the entrypoint, but the build has completed for the image.
if [ -z "$RELEASE_TEST" ]; then
  /entrypoints/build.sh
fi

# If SSH keys directory exists and setup script exists, run it
if [ -d "/etc/ssh-keys-ro" ] && [ -x "/entrypoints/ssh-setup.sh" ]; then
  echo "Running SSH setup script"
  /entrypoints/ssh-setup.sh
  echo "SSH setup script completed"
fi

rm -f tmp/pids/*.pid

if [ -n "$WAIT_FOR_DB" ]; then
  dockerize -wait tcp://${APP_NAME}_db:5432 -timeout 500s
fi

if [ -e /app/development ]; then
  for hook in /app/development/*; do
    [ -x "$hook" ] && $hook
  done
fi

if [ -n "$UPDATE_STATE" ]; then
  ./bin/${APP_NAME} migrate
eval "${app_name_capitalized}_ENV=test ./bin/${APP_NAME} migrate"
fi

mkdir -p /app/data/uploads
mkdir -p /app/data/data_blocks

exec $@
