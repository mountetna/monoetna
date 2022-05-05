include ../make-base/stubs.mk
include ../make-base/docker-compose.mk

app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
app_service_name:=${app_name}_app
fullTag:=$(IMAGES_PREFIX)$(app_name)$(IMAGES_POSTFIX)

run-image-test:: compose-ready
	docker run --rm -e "HOST_DIR=$${PWD}" -v "$${PWD}/tests/outputs:/app/tests/outputs" -v "archimedes-exec:/archimedes-exec" --mount 'type=bind,source=/var/run/docker.sock,target=/var/run/docker.sock,readonly' -e RELEASE_TEST=1 -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} ${fullTag} /entrypoints/development.sh pytest

# Ensures that the etnaagent/archimedes-node image is ready
compose-ready::
	@ IMAGES_PREFIX=etnaagent/ IMAGES_POSTFIX=:latest make -C ../archimedes-node release-build
