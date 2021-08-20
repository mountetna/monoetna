baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

image:
	../../docker/build_image -f Dockerfile
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi
