#!/usr/bin/env bash

set -e

prog=$0
error() {
   echo "Whoops!  Looks like $1:$2 failed."
   exit 1
}
trap 'error "${BASH_SOURCE}" "${LINENO}"' ERR

echoRun() {
  echo "$@"
  $@
}

postToSlack() {
  post-to-slack "deployer on $(hostname)" "watchtower-ping" "$1"
}

findContainers() {
  while read -r line; do
    set -- $line
    containerId="$1"
    labels="$2"

    if echo "$labels" | grep -v "edu.ucsf.mountetna.deployer.disable=true" &>/dev/null; then
      echo "$containerId"
    fi
  done < <(docker ps --format "table {{.ID}}\t{{.Labels}}" $@ | tail -n+2)
}

findImages() {
  docker inspect $(findContainers) | jq '[.[].Image | sub("sha256:"; "")] | unique | .[]' -r
}

findRepoTags() {
  docker inspect $(findImages) | jq '.[] | select(.RepoTags | length > 0) | .RepoTags[]' -r
}

findNonRepoTaggedImages() {
  docker inspect $(findImages) | jq '.[] | select(.RepoTags | length == 0) | .Id' -r
}

pullImages() {
  while read -r line; do
    set -- $line
    if [[ $1 =~ .*/.*:.* ]]; then
      echoRun docker pull $1
    else
      echo "Skipping $1, not remote image"
    fi
  done < <(findRepoTags)
}

linksForContainer() {
  docker inspect  $(docker ps -aq) | jq -r "select(.[].NetworkSettings.Networks | map(.Links[]?) | any(match(\"${1}:\")))"
}

stopNonRepoTaggedImages() {
  while read line; do
    set -- $line
    imageId=$1
    while read line; do
      set -- $line
      containerId=$1
      linkedFailed=0

      postToSlack "Newer version found, attempting to stop $containerId"

      while read line; do
        set -- $line
        linkedContainer=$1
        set +e
        echoRun docker stop $linkedContainer
        if [[ $? -ne 0 ]]; then
          echo "Stopping $linkedContainer failed, skipping stop for $containerId"
          linkedFailed=1
        fi
        set -e
      done < <(linksForContainer $containerId)

      if [[ $linkedFailed -eq 0 ]]; then
        echoRun docker stop $containerId || postToSlack "Could not stop container $containerId, check logs."
      else
        postToSlack "Linked containers failed to be stopped, not stopping $containerId."
      fi
    done < <(findContainers --filter "ancestor=$imageId")
  done < <(findNonRepoTaggedImages)
}

# no while loop here
# allow systemd to kick this process back up automatically, so that we can test it as a single
# iteration.
pullImages
stopNonRepoTaggedImages
docker container prune --force
docker image prune --force
[[ -z "$IS_TEST" ]] && sleep 300 || true
