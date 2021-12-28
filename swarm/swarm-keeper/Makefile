include ../../make-base/stubs.mk

.PHONY: bash
bash:: release-build
	docker run --network host -e SOCKS=$$SOCKS --mount "type=bind,source=$${PWD}/lib/swarm-keeper,target=/usr/lib/swarm-keeper" --rm -it swarm-keeper bash
