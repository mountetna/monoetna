#!/usr/bin/env bash

set -e

prog=$0
error() {
   echo "Whoops!  Looks like $1:$2 failed."
   exit 1
}
trap 'error "${BASH_SOURCE}" "${LINENO}"' ERR

echo "Starting deployer..."

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
  done < <(docker ps --format "table {{.Names}}\t{{.Labels}}" $@ | tail -n+2)
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

containersUsingOutOfDateVolumesFrom() {
  local myMounts="$(docker inspect $1 | jq '[.[].Mounts[].Name]')"
  docker inspect $(docker ps -aq) | jq --argjson mounts "${myMounts}" -r ".[] | select(([.HostConfig.VolumesFrom[]?] | any(. == \"$1\")) and ([.Mounts[].Name] | contains(\$mounts) | not)) | .Name | sub(\"^/\"; \"\")"
}

stopContainer() {
  echoRun docker stop -t ${DOCKER_STOP_TIMEOUT:-10} $1
  eval "${_TEST_HOOK:-true}"
}

# Find containers whose running image does not have a repo tag, and stop them.  Allow systemd to initiate the restart
# so that it is attached to the fd for logging.
# Repo tags move between the images that are most recent for that given tag name.
# Thus, the lack of repo tags indicates being out of date -- the tag has moved to another container.
# The presence of a repo tag indicates that no more recent form of that tag exists on another image.
stopNonRepoTaggedImages() {
  while read line; do
    set -- $line
    imageId=$1
    while read line; do
      set -- $line
      containerId=$1
      linkedFailed=0

      postToSlack "Newer version found, attempting to stop $containerId"
      stopContainer $containerId || postToSlack "Could not stop container $containerId, check logs."
    done < <(findContainers --filter "ancestor=$imageId")
  done < <(findNonRepoTaggedImages)

  while read line; do
    set -- $line
    containerId=$1
    while read line; do
      set -- $line
      dependent=$1

      postToSlack "Container $dependent requires updated $containerId data, restarting."
      stopContainer $dependent || postToSlack "Could not stop container $dependent, check logs."
    done < <(containersUsingOutOfDateVolumesFrom $containerId)
  done < <(findContainers)
}

# no while loop here
# allow systemd to kick this process back up automatically, so that we can test it as a single
# iteration.
pullImages
stopNonRepoTaggedImages
docker container prune --force
docker image prune --force
[[ -z "$IS_TEST" ]] && sleep 300 || true
