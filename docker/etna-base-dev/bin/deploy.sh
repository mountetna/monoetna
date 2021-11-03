#!/bin/bash

set -e

function notifySlackOfVersion() {
  if [ -e /var/run/deployment ]; then
    mkdir -p "/var/run/deployment/notified/${APP_NAME}"
    if ! [ -e "/var/run/deployment/notified/${APP_NAME}/$(cat /built-from-sha)" ]; then
      touch "/var/run/deployment/notified/${APP_NAME}/$(cat /built-from-sha)"
      post-to-slack "${APP_NAME} in ${ENV}" "bioinformatics-ping" "Deployed https://github.com/mountetna/monoetna/commit/$(cat /built-from-sha)"
    fi
  fi
}

function syncAssets() {
  if [ -e /app/public ]; then
    rsync -azvi --delete /app/public/ /sync-assets
  fi
}

function migrate() {
  /app/bin/$APP_NAME migrate
}

function globalMigrate() {
  /app/bin/$APP_NAME global_migrate
}

while [ "$#" -gt 1 ]; do
  if [ "$1" = "--" ]; then
    shift
    break
  fi

  echo "> Running $1"
  $1
  shift
done

exec $@
