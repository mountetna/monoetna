#!/usr/bin/env bash
set -e
#if [ -n "$VERBOSE" ]; then
  set -x
#fi
# Default directories
mkdir -p /app/tmp/pids
mkdir -p /app/public/js
mkdir -p /app/public/css
mkdir -p /app/public/images
mkdir -p /app/log
mkdir -p /app/vendor/bundle
mkdir -p /app/data/
shopt -s globstar
if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle install -j "$(nproc)"
fi
if [ -n "$RUN_NPM_INSTALL" ]; then
  # The images tend to build as root, which for host systems is unsafe,
  # but in containers is fine.
  npm install --unsafe-perm
fi
# Prepare the httpd conf hooks.  These are shared from the app images to the app_fe image.
# This allows individual apps to package their frontend configuration but run them
# as a separate docker container.
mkdir -p /usr/opt/httpd.conf.d
echo "@app_name = '${APP_NAME}'" > /usr/opt/vars.rb
rm -rf /usr/opt/httpd.conf.d/*.include
if [ -e /app/build ]; then
  for hook in /app/build/*; do
    if stat -c  %A $hook | grep x &>/dev/null; then
      $hook
    fi
  done
fi