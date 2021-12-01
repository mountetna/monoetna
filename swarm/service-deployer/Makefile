include ../../make-base/stubs.mk

.PHONY: bash
bash:: image
	docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock:ro $$fullTag bash
