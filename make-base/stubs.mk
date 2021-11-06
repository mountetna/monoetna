include $(shell dirname $(lastword $(MAKEFILE_LIST)))/utils.mk

# Set this up as the default context for any Makefile.
# Every project's default name matches the directory it exists in,
# as per utils.mk's find_project_containing logic.
app_name:=$(shell basename $$(pwd))
fullTag=$(IMAGES_PREFIX)$(app_name)$(IMAGES_POSTFIX)

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

.PHONY: restart
restart:: ## Restarts all containers for the given project(s)
	@ true

.PHONY: release
release:: release-build release-test ## Builds static docker images staged for release, runs tests against them, and pushes them to dockerhub (requires PUSH_IMAGES=1)
	@ true

.PHONY: release-build
release-build:: config-ready ## Step that prepares a release build if necessary
	@ true

.PHONY: release-test ## Step that prepares a release test if necessary
release-test:: config-ready release-build
	@ true

.PHONY: config-ready
config-ready:: ## Setup step that ensures that configuration files necessary for preparation of the development environment are ready.
	@ true

.PHONY: update
update:: update-ready ## Step to force update of database and development dependencies.
	@ true

.PHONY: update-ready
update-ready::
	@ true

.PHONY: compose-ready
compose-ready:: ## Setup step that ensures that the images are ready to run the development docker-compose.yml file
	@ true

.PHONY: _test-find-all
_test-find-all:
	@ echo $(call map,dirname,$(call find_projects_containing,Makefile))

ifneq ($(wildcard ./Dockerfile),)
$(eval $(call image_target_here))
endif
