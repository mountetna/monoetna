#!/usr/bin/env bash

set -e
set -o pipefail

error() {
   echo "Whoops!  Looks testing failed at $1:$2"
   exit 1
}
trap 'error "${BASH_SOURCE}" "${LINENO}"' ERR

assertRunning() {
  if docker ps --filter name=$1 | grep $1 &>/dev/null; then
    return 0
  fi

  return 1
}

assertNotRunning() {
  if docker ps --filter name=$1 | grep -v $1 &>/dev/null; then
    return 0
  fi

  return 1
}

docker kill test_a || true
docker kill test_b || true
docker kill test_c || true
docker kill test_d || true

docker image rm image_a:latest || true
docker image rm image_b:latest || true

docker build --rm --force-rm --tag image_a:latest --build-arg FILE=a -f TestDockerFile /root
docker build --rm --force-rm --tag image_b:latest --build-arg FILE=b -f TestDockerFile /root

docker run -d --name test_a --rm image_a sleep 100
docker run -d --name test_b --rm image_a sleep 100
docker run -d --name test_c --rm image_b sleep 100
docker run -d --name test_d --rm image_b sleep 100

assertRunning test_a
assertRunning test_b
assertRunning test_c
assertRunning test_d

export IS_TEST=1
deployer

# Does not stop any containers with up to date images.
assertRunning test_a
assertRunning test_b
assertRunning test_c
assertRunning test_d

docker build --rm --force-rm --tag image_a:latest  --build-arg FILE=c -f TestDockerFile /root

deployer

# Does not stop any containers with up to date images.
assertNotRunning test_a
assertNotRunning test_b
assertRunning test_c
assertRunning test_d

docker run -d --name test_a --rm --link test_c:somename image_a sleep 100

deployer

assertRunning test_a
assertNotRunning test_b
assertRunning test_c
assertRunning test_d

docker build --rm --force-rm --tag image_b:latest  --build-arg FILE=d -f TestDockerFile /root

deployer

assertNotRunning test_a
assertNotRunning test_b
assertNotRunning test_c
assertNotRunning test_d
