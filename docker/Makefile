export COMPOSE_PROJECT_NAME=monoetna
projects := $(shell ls ../*/Makefile | grep -v docker | xargs -n 1 dirname | xargs -n 1 basename)
compose_ymls     := $(shell ls ../*/docker-compose.yml | grep -v 'docker/docker-compose.yml')
images := $(shell echo */)

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

docker-compose.yml: $(compose_ymls) ./compose $(wildcard ./*.shared.yml) ./default_compose
	@ for project in $(projects); do cd ../$$project && make docker-compose.yml; done
	@ echo 'Rebuilding docker-compose.yml from projects $(compose_ymls)'
	@ ./compose $(compose_ymls) > .tmp
	@ mv .tmp docker-compose.yml

.PHONY: config-ready
config-ready: docker-compose.yml ## Prepares the top level docker-compose.yml config.  Run this if you recently changed docker-compose.yml files.

.PHONY: docker-ready
docker-ready: config-ready build

.PHONY: build
build: ## Builds all development related docker images
	@ for image in $(images); do make -s -C $$image || exit; done

.PHONY: update
update: docker-ready ## Runs all projects' updates, running bundle and npm install
	@ set -e && for project in $(projects); do make -C ../$$project update; done

.PHONY: up
up: docker-ready ## Starts up all containers of this project in the background
	@ for project in $(projects); do cd ../$$project && make config.yml || true; done
	@ docker-compose up -d

.PHONY: down
down: docker-compose.yml ## Ends all projects' processes
	@ docker-compose down --remove-orphans

.PHONY: ps
ps: ## Shows ps of all projects' containers
	@ docker-compose ps

.PHONY: logs
logs: ## Shows logs of all running projects' containers
	@ docker-compose logs -f

.PHONY: logs-recent
logs-recent: ## For CI
	@ docker-compose logs

.PHONY: clean
clean: docker-ready
	@ docker-compose down -v --remove-orphans
	@ docker-compose create
