#!/bin/bash
set -e

if [ -z "$SKIP_RUBY_SETUP" ]; then
  bundle check || bundle install -j "$(nproc)"
  rm -f tmp/pids/*.pid

  mkdir -p tmp/pids

  if [ -z "$SKIP_DB_WAIT" ]; then
    dockerize -wait tcp://magma_db:5432 -timeout 60s

    for project in $(findEnvConfig | grep ':project_path:' | sed -e 's/.*:project_path:\s*//g'); do
      echo "Initializing $project..."
      bin/magma create_db "$(basename $project)"
    done

    bin/magma global_migrate
    bin/magma migrate
  fi
fi

exec "$@"
