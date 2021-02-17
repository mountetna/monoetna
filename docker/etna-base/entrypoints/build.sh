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
#  if [ -e ../etna ]; then
#    bundle config set disable_local_branch_check true
#    bundle config set local.etna /etna
#  fi
  # Clear any existing checkouts that occurs in development volumes.
  rm -rf /app/vendor/bundle/*/ruby/*/bundler/gems/monoetna-*
  bundle install -j "$(nproc)"

  # This jank due to the wacky way bundle enforces git repo behavior means we have to hijack the production
  # checkout.  When we can have deploys use the gem inside the container _rather than_ the external bundled one,
  # this will become simpler.
  if [ -e ../etna ] && grep 'monoetna.git' /app/Gemfile 1>/dev/null; then
    mkdir -p /tmp/
    rm -rf /tmp/.git
    mv /app/vendor/bundle/*/ruby/*/bundler/gems/monoetna-*/.git /tmp/.git
    rm -rf /app/vendor/bundle/*/ruby/*/bundler/gems/monoetna-*/*
    cp -r /etna/* /app/vendor/bundle/*/ruby/*/bundler/gems/monoetna-*/
    rm -rf /app/vendor/bundle/*/ruby/*/bundler/gems/monoetna-*/packages
    mv /tmp/.git /app/vendor/bundle/*/ruby/*/bundler/gems/monoetna-*/
    cd /app
    bundle install
    rm /app/vendor/bundle/*/bin/etna
  fi
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

if ! [ -e config.yml ] && [ -e config.yml.template ]; then
  cp config.yml.template config.yml
fi

if [ -n "$RUN_NPM_INSTALL" ]; then
  if [ -e ../etna ]; then npm link ../etna/packages/etna-js; fi
fi

if [ -n "$RELEASE_TEST" ]; then
  # TODO: Find a better test harness for this.  In release tests, let's also ensure that in the resulting
  # setup, the app_fe comes up successfully.
  dockerize -wait tcp://${APP_NAME}_app_fe:80 -timeout 300s
fi

echo 'for file in /app/*.completion; do source $file || true; done' >> /root/.bashrc
echo 'export PATH="/app/bin:$PATH"' >> /root/.bashrc
# Allow other users to use the root bash setup
chmod 744 -R /root/
