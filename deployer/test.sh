#!/usr/bin/env bash

set -e
set -o pipefail

export DOCKER_STOP_TIMEOUT=1

error() {
   echo "Whoops!  Looks testing failed at $1:$2"
   exit 1
}
trap 'error "${BASH_SOURCE}" "${LINENO}"' ERR

red="\e[0;91m"
green="\e[0;92m"
bold="\e[1m"
reset="\e[0m"

assert() {
  if $@; then
    echo -e ${green}${bold}$@${reset}
    return 0
  fi

  echo -e ${red}${bold}not $@${reset}
  return 1;
}

assertNot() {
  if ! $@; then
    echo -e ${green}${bold}not $@${reset}
    return 0
  fi

  echo -e ${red}${bold}$@${reset}
  return 1;
}

running() {
  if docker ps --filter name=$1 | grep $1 &>/dev/null; then
    return 0
  fi

  return 1
}

imagesOfContainer() {
  docker inspect $1 | jq '[.[].Image | sub("sha256:"; "")] | unique | .[]' -r
}

findNonRepoTaggedImages() {
  docker inspect $(imagesOfContainer $1) | jq '.[] | select(.RepoTags | length == 0) | .Id' -r
}

usingLatest() {
  if [[ "$(findNonRepoTaggedImages $1 2>/dev/null | wc -l)" -eq 0 ]]; then
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

assert usingLatest test_a
assert usingLatest test_b
assert usingLatest test_c
assert usingLatest test_d

assert running test_a
assert running test_b
assert running test_c
assert running test_d

export IS_TEST=1
deployer

# Does not stop any containers with up to date images.
assert running test_a
assert running test_b
assert running test_c
assert running test_d

docker build --rm --force-rm --tag image_a:latest  --build-arg FILE=c -f TestDockerFile /root

assertNot usingLatest test_a
assertNot usingLatest test_b
assert usingLatest test_c
assert usingLatest test_d

deployer

# Does not stop any containers with up to date images.
assertNot running test_a
assertNot running test_b
assert running test_c
assert running test_d

RESTART_TEST_A="docker run -d --name test_a --rm --volumes-from test_c image_a sleep 100"
$RESTART_TEST_A

deployer

assert running test_a
assertNot running test_b
assert running test_c
assert running test_d

docker build --rm --force-rm --tag image_b:latest  --build-arg FILE=d -f TestDockerFile /root

#_TEST_HOOK="$RESTART_TEST_A 2>/dev/null" deployer
#
#assert running test_a
#assertNot running test_b
#assertNot running test_c
#assertNot running test_d
#
#assert usingLatest test_a
