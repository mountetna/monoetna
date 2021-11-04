include $(shell dirname $(lastword $(MAKEFILE_LIST)))/utils.mk
include $(shell dirname $(lastword $(MAKEFILE_LIST)))/build.mk

# Set this up as the default context for any Makefile.
# Every project's default name matches the directory it exists in,
# as per utils.mk's find_project_containing logic.
app_name:=$(shell basename $$(pwd))
fullTag=$(IMAGES_PREFIX)$(app_name)$(IMAGES_POSTFIX)
fullFeTag=$(IMAGES_PREFIX)$(app_name)_app_fe$(IMAGES_POSTFIX)

# a phony target that always forces build of targets that depend on it.
.PHONY: always

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | uniq | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

# This file is intended to be scoped to a single project root at a time.
# Each command is imperative, it invokes a process and can be extended to
# define a process for a project scope.
# Consider these as entry points for make -C <project> <thing-to-do>

.PHONY: up
up:: ## Starts up the docker containers associated with the given project(s) in the background
	@ true

.PHONY: down
down:: ## Ends the background docker containers associated with the given project(s)
	@ true

.PHONY: ps
ps:: ## Shows status of the containers running with the given project(s)
	@ true

.PHONY: logs
logs:: ## Shows logs of the running containers with the given project(s)
	@ true

.PHONY: bash
bash:: ## Starts a bash shell in an app environment for the given project
	@ true

.PHONY: psql
psql:: ## Starts a psql shell in an app environment for the given project
	@ true

.PHONY: migrate
migrate:: ## Runs migrations in an app environment for the given project
	@ true

.PHONY: restart
restart:: ## Restarts all containers for the given project(s)
	@ true

.PHONY: release
release:: release-build release-test ## Builds static docker images staged for release, runs tests against them, and pushes them to dockerhub (requires PUSH_IMAGES=1)
	if docker inspect --type=image $(fullTag) &>/dev/null; then docker push $(fullTag); fi
	if docker inspect --type=image $(fullFeTag) &>/dev/null; then docker push $(fullFeTag); fi

.PHONY: release-build
release-build:: config-ready ## Step that just builds docker images staged for release
	@ true

.PHONY: release-test ## Step that just runs the tests for currently built release images.
release-test::
	@ true

config-ready:: ## Setup step that ensures that configuration files necessary for preparation of the development environment are ready.
	@ true

update:: ## Step to force update of database and development dependencies.
	@ true

compose-ready:: ## Setup step that ensures that the images are ready to run the development docker-compose.yml file
	@ true

.PHONY: _test-find-all
_test-find-all:
	@ echo $(call map,dirname,$(call find_projects_containing,Makefile))
