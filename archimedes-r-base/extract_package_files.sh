#!/usr/bin/env bash

set -e
set -x

base='archimedes-r-base'
prefix='etnaagent/'
# This should only ever be run on the :master image.
# By the time anything is pushed to :production, this should already have been done!
tag=':master'

target="${prefix}${base}${tag}"

docker pull "$target"

container_id=$(docker create "$target")
docker cp "$container_id:/app/DESCRIPTION" .
docker cp "$container_id:/app/renv.lock" .
docker cp "$container_id:/app/.Rprofile" .
docker rm "$container_id"