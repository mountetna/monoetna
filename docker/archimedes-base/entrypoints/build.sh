#!/usr/bin/env bash

set -e

#if [ -n "$VERBOSE" ]; then
  set -x
#fi

shopt -s globstar

if [ -z "$SKIP_PYTHON_SETUP" ]; then
  poetry install
fi

if [ -e /app/build ]; then
  for hook in /app/build/*; do
    [ -x "$hook" ] && $hook
  done
fi

if ! [ -e config.yml ] && [ -e config.yml.template ]; then
  cp config.yml.template config.yml
fi

poetry completions bash > /app/poetry.completion
touch /root/.bashrc
echo 'for file in /app/*.completion; do source $file || true; done' >> /root/.bashrc
echo 'export PATH="/app/bin:$PATH"' >> /root/.bashrc
# Allow other users to use the root bash setup
chmod -R 744 /root/
