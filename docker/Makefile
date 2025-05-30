include ../make-base/utils.mk
export COMPOSE_PROJECT_NAME=monoetna
compose_ymls:=$(filter-out %/docker/docker-compose.yml,$(call find_projects_containing,docker-compose.yml))
compose_projects:=$(call map,dirname,$(compose_ymls))

.PHONY: help
help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help

docker-compose.yml: $(compose_ymls) ./build_support/compose $(wildcard ./*.shared.yml)
	@ echo 'Rebuilding docker-compose.yml from projects $(compose_ymls)'
	@ ./build_support/compose $(compose_ymls) > .tmp
	@ mv .tmp docker-compose.yml

.PHONY: config-ready
config-ready: docker-compose.yml ## Prepares the top level docker-compose.yml config.  Run this if you recently changed docker-compose.yml files.

.PHONY: compose-ready
compose-ready: config-ready

.PHONY: update
update:

.PHONY: update-all
update-all: compose-ready ## Runs all projects' updates, running bundle and npm install
	@ set -e && for project in $(compose_projects); do make -C $$project update; done

.PHONY: up
up: compose-ready ## Starts up all containers of this project in the background
	@ set -e; for project in $(compose_projects); do make -C $$project compose-ready; done
	@ docker compose up -d

.PHONY: down
down: docker-compose.yml ## Ends all projects' processes
	@ docker compose down --remove-orphans

.PHONY: ps
ps: ## Shows ps of all projects' containers
	@ docker compose ps

.PHONY: logs
logs: ## Shows logs of all running projects' containers
	@ docker compose logs -f

.PHONY: logs-recent
logs-recent: ## For CI
	@ docker compose logs
