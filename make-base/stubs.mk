# Set this up as the default context for any Makefile.
# Every project's default name matches the directory it exists in,
# as per build.mk's find_project_containing logic.
app_name:=$(shell basename $$(pwd))

# Function that find project directories by a file that must exist inside of it
define find_projects_containing
$(sort \
	$(call map,readlink, \
		$(filter-out $(shell dirname $(firstword $(wildcard $(addsuffix /monoetna,.. ../.. ../../..))))/%, \
			$(wildcard \
				$(addsuffix /*/$(1),. .. ../.. ../../..) \
				$(addsuffix /docker/*/$(1),. .. ../.. ../../..) \
				$(addsuffix /swarm/*/$(1),. .. ../.. ../../..) \
			) \
		) \
	) \
)
endef

define find_project_containing
$(filter %/$(1)/$(2),$(call find_projects_containing,$(2)))
endef

map = $(foreach a,$(2),$(call $(1),$(a)))
readlink = $(shell readlink -m $(1))
dirname = $(shell basename $(shell dirname $(1)))
mkfile_path=$(abspath $(lastword $(MAKEFILE_LIST)))

make_base:=$(call find_project_containing,make-base,build.mk.template)

# a phony target that always forces build of targets that depend on it.
.PHONY: always
always:

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

.PHONY: build
build:: ## Forces a rebuild of project development dockerfiles
	@ true

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
release:: ## Builds static docker images staged for release, runs tests against them, and pushes them to dockerhub (requires PUSH_IMAGES=1)
	@ true

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

.PHONY: docker_buildable
docker_buildable: Dockerfile build.mk

Dockerfile:
	@ echo '# Check other project Dockerfiles out for how to set this up' > Dockerfile

# Creates the project specific global namespaced build file
# Always re-evaluate this to ensure up to date
build.mk: always
	$(make_base)/../bin/bash-template $(make_base)/build.mk.template > build.mk

.PHONY: _test-find-all
_test-find-all:
	@ echo $(call map,dirname,$(call find_projects_containing,Makefile))
