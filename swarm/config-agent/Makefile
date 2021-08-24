baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

.PHONY: image
image:
	../../docker/build_image -f Dockerfile
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi

.PHONY: test
test: image
	docker run --rm -v /var/run/docker.sock:/var/run/docker.sock:ro -v $$PWD/test:/script/ config-agent