include ../make-base/stubs.mk

app_name=archimedes
include ../make-base/docker-compose.mk

app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
app_service_name:=${app_name}_app
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

release-test:: docker-ready
	docker run --rm -e "HOST_DIR=$${PWD}" -v "${PWD}/tests/outputs:/app/tests/outputs" -v "archimedes-exec:/archimedes-exec" --mount 'type=bind,source=/var/run/docker.sock,target=/var/run/docker.sock,readonly' -e RELEASE_TEST=1 -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} ${fullTag} /entrypoints/development.sh pytest