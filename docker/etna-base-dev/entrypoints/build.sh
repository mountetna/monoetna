#!/usr/bin/env bash

set -e
set -x

# Default directories
mkdir -p /app/tmp/pids
mkdir -p /app/public/js
mkdir -p /app/public/css
mkdir -p /app/public/images
mkdir -p /app/log
mkdir -p /app/vendor/bundle
mkdir -p /app/data/

shopt -s globstar

if [ -n "$FULL_INSTALL" ]; then
  bundle install -j "$(nproc)"
fi

if [ -n "$FULL_INSTALL" ]; then
  # The images tend to build as root, which for host systems is unsafe,
  # but in containers is fine.
  npm install --unsafe-perm
fi

# Prepare the puma.sh into the /app/bin directory.
ln -sf /entrypoints/puma.sh /app/bin/
chmod +x /app/bin/puma.sh

if [ -e /app/build ]; then
  for hook in /app/build/*; do
    if stat -c  %A $hook | grep x &>/dev/null; then
      $hook
    fi
  done
fi

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
