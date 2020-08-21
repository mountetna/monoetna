CERTS ?= ../../development-certs/certs
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

image:
	if [ -n "$$PULL_IMAGES" ]; then docker pull $(fullTag); fi
	../build_image Dockerfile $(CERTS) ../../etna
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi
