CERTS ?= ../../development-certs/certs
DOCKER ?= ../../docker
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

image:
	../build_image Dockerfile $(CERTS) $(DOCKER)
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi
