help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | uniq | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

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
release:: ## Builds static docker images staged for release, runs tests against them, and pushes them to dockerhub (requires PUSH=1)
	@ true

.PHONY: release-build
release-build:: config-ready
	@ true

.PHONY: release-test
release-test::
	@ true

.PHONY: irb
irb:: ## Starts up an irb session in the context of the given project
	@ true

config-ready::
	@ true

update::
	@ true

docker-ready::
	@ true
