uid=$(shell id -u ${USER})
gid=$(shell getent group docker | cut -d':' -f3)
volumes:=-v $$PWD:/data -v /var/run/docker.sock:/var/run/docker.sock -v /etc/passwd:/etc/passwd:ro -v /etc/group:/etc/group:ro -v /var/mountetna:/var/mountetna
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

environments/local.json:
	cp environments/local.json.template environments/local.json

.PHONY: create
create: images environments/local.json ## Creates and converges a running centos image with the kitchen.yml settings
	$(chefdk) kitchen converge

.PHONY: ssh
ssh: images ## Connects via ssh to the running centos image created last
	$(chefdk) ssh -i .kitchen/docker_id_rsa kitchen@$$(docker inspect centos-7 | grep IPAddress | tail -n1 | cut -d'"' -f4) -p 22

.PHONY: destroy
destroy: images ## destroys the running test centos deployment.  Note, `docker container prune` may be necessary if there is a conflict.
	$(chefdk) kitchen destroy
	sudo rm -rf /var/mountetna/*

.PHONY: bash
bash: images ## Opens a bash environment inside the chefdk tool
	$(chefdk)
