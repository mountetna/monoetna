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
mkdir -p /tmp/metrics.prom/

# Often, etna is built as root, but the processes are not run as root.
# Certain, and only certain, directories should be writable.
# Especially with source code in script-based servers, we don't want to
# generally allow write-ability into the source code by any logic.
chmod -R 777 /app/tmp
chmod -R 777 /app/log
chmod -R 777 /tmp/metrics.prom/

shopt -s globstar

if [ -n "$FULL_BUILD" ]; then
  bundle install -j "$(nproc)"
fi

if [ -n "$FULL_BUILD" ]; then
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
