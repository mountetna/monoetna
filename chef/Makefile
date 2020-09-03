uid=$(shell id -u ${USER})
gid=$(shell getent group docker | cut -d':' -f3)
volumes:=-v $$PWD:/data -v /var/run/docker.sock:/var/run/docker.sock -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro
chefdk:=docker run --rm -ti -u $(uid):$(gid) $(volumes) test-chefdk

# ssh -v -i .kitchen/docker_id_rsa kitchen@172.17.0.3 -p 22
help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

.PHONY: images
images:
	make -C ../docker build

.PHONY: test
test: images ## Initializes a docker using centos and the chef recipes referenced in kitchen.yml, testing the deployment.
	$(chefdk) kitchen test

.PHONY: create
create: images ## Creates a test centos deployment using docker
	$(chefdk) kitchen create

.PHONY: destroy
destroy: images ## destroys the running test centos deployment
	$(chefdk) kitchen destroy

.PHONY: bash
bash: images ## Opens a bash environment inside the chefdk tool
	$(chefdk)
