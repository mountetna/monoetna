#!/bin/bash

set -e

function notifySlackOfVersion() {
  if [ -e /var/run/deployment ]; then
    mkdir -p "/var/run/deployment/notified/${APP_NAME}"
    if ! [ -e "/var/run/deployment/notified/${APP_NAME}/$(cat /built-from-sha)" ]; then
      touch "/var/run/deployment/notified/${APP_NAME}/$(cat /built-from-sha)"
      post-to-slack "${APP_NAME} in ${ENV}" "bioinformatics-ping" "Deployed https://github.com/mountetna/monoetna/commit/$(cat /built-from-sha)"

      if [[ $APP_NAME == "metis" && $ENV == "production" ]]; then
        post-to-slack "${APP_NAME} in ${ENV}" "bioinformatics-ping" "<!channel> ${APP_NAME} restarted in ${ENV} -- please check your upload / download jobs!"
      fi
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
