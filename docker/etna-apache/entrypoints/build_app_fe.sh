#!/usr/bin/env bash

set -e
set -x

shopt -s globstar

PROXY_APP_NAME=
if [[ "${APP_NAME}" =~ (.*)_app_fe$ ]]; then
  PROXY_APP_NAME=${BASH_REMATCH[1]}
else
  exit 1
fi


# Prepare the httpd conf hooks.  These are shared from the app images to the app_fe image.
# This allows individual apps to package their frontend configuration but run them
# as a separate docker container.
mkdir -p /usr/opt/httpd.conf.d
# Clear any existing as part of some previous build.
rm -rf /usr/opt/httpd.conf.d/*.include
sed main.conf.template -e "s/-app_name-/$PROXY_APP_NAME/g" > /usr/opt/httpd.conf.d/main.conf

if [ -e "/app/build" ]; then
  for hook in /app/build/*; do
    if stat -c  %A $hook | grep x &>/dev/null; then
      $hook
    fi
  done
fi
