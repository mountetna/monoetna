include ../../make-base/stubs.mk

baseTag:=$(shell basename "$$(pwd)")

.PHONY: bash
bash:: release-build
	docker run --rm -it -v ${PWD}:/app $(baseTag) bash

.PHONY: metis
metis:: release-build
	docker run --rm -it -v ${PWD}/..:/swarm $(baseTag) npx ts-node apply-yaml.ts /swarm/metis.yml
