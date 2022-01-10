include ../../make-base/stubs.mk

fullTag:=$(IMAGES_PREFIX)$(app_name)$(IMAGES_POSTFIX)

.PHONY: bash
bash:: release-build
	docker run --network host -e SOCKS=$$SOCKS --mount "type=bind,source=$${PWD}/lib/swarm-keeper,target=/usr/lib/swarm-keeper" --rm -it swarm-keeper bash

run-image-test::
	docker run --rm -e ENDPOINT=1 -e SOCKS=127.0.0.1:9090 -e RELEASE_TEST=1 -e IS_CI=$${IS_CI} ${fullTag} /usr/lib/swarm-keeper-tests
