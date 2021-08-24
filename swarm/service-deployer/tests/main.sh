#!/usr/bin/env bash
set -euo pipefail
source /bach.sh

test-watch-and-update-service() {
  @mock docker service ls --format "{{.Name}}" --filter label=autoupdate=true \
    === @stdout service-c \
    service-b \
    service-a

  @mock docker service inspect service-a --format "{{.Spec.TaskTemplate.ContainerSpec.Image}}" \
    === @stdout coolguy/image:latest@sha256:a
  @mock docker service inspect service-b --format "{{.Spec.TaskTemplate.ContainerSpec.Image}}" \
    === @stdout coolguy/image2:latest@sha256:b
  @mock docker service inspect service-c --format "{{.Spec.TaskTemplate.ContainerSpec.Image}}" \
    === @stdout coolguy/image:staging@sha256:c

  @mock docker image inspect sha256:a \
    === @stdout <<-EOF
[{
        "RepoTags": []
}]
EOF

  @mock docker image inspect sha256:b \
    === @stdout <<-EOF
[{
        "RepoTags": ["not-it"]
}]
EOF

  @mock docker image inspect sha256:c \
    === @stdout <<-EOF
[{
        "RepoTags": ["coolguy/image:staging"]
}]
EOF

  @@mock

  @ignore echo
  @run /subject/watch-and-update-service-images
}

grep() {
  @real grep
}

test-watch-and-update-service-assert() {
  echo "Starting deployer..."
  echo ""
}