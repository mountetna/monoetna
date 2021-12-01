export COMPOSE_PROJECT_NAME=monoetna
images := $(shell echo */)

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

.PHONY: build
build: ## Builds all swarm related docker images
	@ for image in $(images); do make -s -C $$image || exit; done
