include ../make-base/stubs.mk

EXTRA_DOCKER_ARGS:=-v /var/run/docker.sock:/var/run/docker.sock:ro $(EXTRA_DOCKER_ARGS)
app_name=vulcan
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk
include ../make-base/node.mk

# Ensures that the etnaagent/archimedes image is ready
docker-ready::
	@ IMAGES_PREFIX=etnaagent/ make -C ../archimedes release-build
