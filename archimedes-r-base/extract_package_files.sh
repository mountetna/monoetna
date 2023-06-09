#!/usr/bin/env bash

set -e
set -x

base='archimedes-r-base'
prefix='' # 'etnaagent/' if should pull from remote
# This should only ever be run on the :latest or :master image.
# By the time anything is pushed to :production, this should already have been done!
tag=':latest'

target="${prefix}${base}${tag}"

container_id=$(docker create "$target")
docker cp "$container_id:/app/DESCRIPTION" .
docker cp "$container_id:/app/renv.lock" .
docker cp "$container_id:/app/.Rprofile" .
docker rm "$container_id"