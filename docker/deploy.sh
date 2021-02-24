#!/bin/bash

set -e

function notifySlackOfVersion() {
  if [ -e /var/run/deployment ]; then
    if ! [ -e /var/run/deployment/notified/$(cat /built-from-sha) ]; then
      touch /var/run/deployment/notified/$(cat /built-from-sha)
      post-to-slack "${APP_NAME} in ${ENV}" "bioinformatics-ping" "Deploying https://github.com/mountetna/monoetna/compare/${ENV}...$(cat /built-from-sha)"
    fi
  fi
}

function migrate() {
  /app/bin/$APP_NAME migrate
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
