include ../make-base/stubs.mk

app_name=prometheus
include ../make-base/docker-compose.mk
EXTRA_DOCKER_ARGS:=-v /var/run/docker.sock:/var/run/docker.sock:ro

containerSh:=sh
