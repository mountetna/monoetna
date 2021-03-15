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

# Not all apache directives support environment variable interpolation, so the best way to hook in
# certain per app variables is via a template.
erb -r /usr/opt/vars.rb -- /opt/fe.conf.erb > /usr/opt/httpd.conf.d/main.conf

if ! [ -e config.yml ] && [ -e config.yml.template ]; then
  cp config.yml.template config.yml
fi

if [ -n "$RELEASE_TEST" ]; then
  # TODO: Find a better test harness for this.  In release tests, let's also ensure that in the resulting
  # setup, the app_fe comes up successfully.
  dockerize -wait tcp://${APP_NAME}_app_fe:80 -timeout 300s
fi

echo 'for file in /app/*.completion; do source $file || true; done' >> /root/.bashrc
echo 'export PATH="/app/bin:$PATH"' >> /root/.bashrc
# Allow other users to use the root bash setup

chmod -R 744 /root
