CERTS ?= ../../development-certs/certs
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

image:
	../build_image Dockerfile $(CERTS)
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi
