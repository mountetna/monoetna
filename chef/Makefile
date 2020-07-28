uid=$(shell id -u ${USER})
gid=$(shell id -g ${USER})

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

markers:
				@ mkdir markers

markers/chef-dk: markers Dockerfile.chefdk
				docker build --tag chefdk-docker . --file Dockerfile.chefdk
				@ touch markers/chef-dk

markers/centos: markers Dockerfile.centos
				docker build --tag centos-docker . --file Dockerfile.centos
				@ touch markers/centos

.PHONY: test
test: markers/chef-dk markers/centos ## Starts up all containers of this project in the background
	@ docker run --rm -ti -v $$PWD:/data -v /var/run/docker.sock:/var/run/docker.sock chefdk-docker kitchen test

.PHONY: create
create: markers/chef-dk markers/centos ## Starts up all containers of this project in the background
	docker run --rm -ti -v $$PWD:/data -v /var/run/docker.sock:/var/run/docker.sock chefdk-docker kitchen create

.PHONY: destroy
destroy: markers/chef-dk markers/centos ## Starts up all containers of this project in the background
	docker run --rm -ti -v $$PWD:/data -v /var/run/docker.sock:/var/run/docker.sock chefdk-docker kitchen destroy

.PHONY: bash
bash: markers/chef-dk markers/centos ## Starts up all containers of this project in the background
	docker run --rm -ti -u $(uid):$(gid) -v $$PWD:/data -v /var/run/docker.sock:/var/run/docker.sock chefdk-docker
