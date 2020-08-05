#!/usr/bin/env sh

set -e

if [ -n "$VERBOSE" ]; then
  set -x
fi

export PATH="/app/node_modules/.bin:/app/vendor/bundle/$RUBY_VERSION/bin:$PATH"
# Default directories
mkdir -p /app/tmp/pids
mkdir -p /app/public/js
mkdir -p /app/public/css
mkdir -p /app/public/images
mkdir -p /app/log
mkdir -p /app/vendor/bundle
mkdir -p /app/data/


if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
fi


if [ -n "$RUN_NPM_INSTALL" ]; then
  npm install --unsafe-perm
fi

# Prepare the httpd conf hooks
mkdir -p /usr/opt/httpd.conf.d
echo "@app_name = '${APP_NAME}'" > /usr/opt/vars.rb
rm -rf /usr/opt/httpd.conf.d/*.include

# NPMRC helps with npm link.
cp /opt/npmrc /app/.npmrc

if [ -e /app/build ]; then
  for hook in /app/build/*; do
    [ -x "$hook" ] && $hook
  done
fi

erb -r /usr/opt/vars.rb -- /opt/fe.conf.erb > /usr/opt/httpd.conf.d/main.conf

