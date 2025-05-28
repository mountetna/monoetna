include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
include ../make-base/etna-ruby.mk
include ../make-base/node.mk
EXTRA_DOCKER_ARGS:=-v /var/run/docker.sock:/var/run/docker.sock:ro -v $${PWD}/../vulcan/spec/data:/app/spec/data -e HOST_DIR_MAP=/app=$${PWD}/../vulcan $(EXTRA_DOCKER_ARGS)
