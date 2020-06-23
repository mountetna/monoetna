export COMPOSE_PROJECT_NAME=monoetna
projects := $(shell ls ../*/Makefile | grep -v docker | xargs -n 1 dirname | xargs -n 1 basename)
compose_ymls     := $(shell ls ../*/docker-compose.yml | grep -v 'docker/docker-compose.yml')

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

docker-compose.yml: $(compose_ymls) ./compose ./empty-docker-compose.yml
				@ echo 'Rebuilding docker-compose.yml from projects $(compose_ymls)'
				@ ./compose $(compose_ymls) > .tmp
				@ mv .tmp docker-compose.yml

.docker-build-mark: $(wildcard ../*/docker/**/*) docker-compose.yml
				$(MAKE) build
				@ touch .docker-build-mark

.PHONY: build
build: docker-compose.yml ## Forces a rebuild of all projects' development dockerfiles
				@ echo $(projects)
				@ set -e; for project in $(projects); do cd ../$$project && make build; done
				@ touch .docker-build-mark

.PHONY: up
up: docker-compose.yml .docker-build-mark ## Starts up all containers of this project in the background
				@ for project in $(projects); do cd ../$$project && make config.yml || true; done
				@ docker-compose up -d

.PHONY: down
down: docker-compose.yml ## Ends all projects' processes
				@ docker-compose down --remove-orphans

.PHONY: ps
ps: docker-compose.yml ## Shows ps of all projects' containers
				@ docker-compose ps

.PHONY: logs
logs: docker-compose.yml ## Shows logs of all running projects' containers
				@ docker-compose logs
