baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

.PHONY: image
image:
	../../docker/build_image -f Dockerfile ../../docker/etna-base-dev
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi

.PHONY: bash
bash: image
	docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock:ro $$fullTag bash
