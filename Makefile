export COMPOSE_PROJECT_NAME=monoetna
projects := $(shell ls ./*/Makefile | grep -v docker | xargs -n 1 dirname | xargs -n 1 basename)
compose_ymls     := $(shell ls ./*/docker-compose.yml)

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

.PHONY: build
build: ## Forces a rebuild of all projects' development dockerfiles
				@ make -C docker build

.PHONY: up
up: ## Starts up all containers of this project in the background
				@ make -C docker up

.PHONY: down
down: ## Ends all projects' processes
				@ make -C docker down

.PHONY: ps
ps: ## Shows ps of all projects' containers
				@ make -C docker ps

.PHONY: logs
logs: ## Shows logs of all running projects' containers
				@ make -C docker logs

.PHONY: bash
bash: ## Starts a bash shell in an app environment
				@ echo Run this within a specific app context, ie: make -C metis bash

.PHONY: psql
psql: ## Starts a psql shell in an app environment
				@ echo Run this within a specific app context, ie: make -C janus psql

.PHONY: migrate
migrate: ## Runs migrations in a specific app context
				@ echo Run this within a specific app context, ie: make -C janus migrate

.PHONY: restart
restart: ## Restarts running services for a project
				@ echo Run this within a specific app context, ie: make -C metis restart

.PHONY: test
test: ## Runs all projects' tests
				@ set -e && for project in $(projects); do make -C $$project test; done

.PHONY: setup-links
setup-links: ## Sets up local development links for npm (and one day ruby) packages to use local copies of files.
				@ for project in $(projects); do make -C $$project setup-links; done
