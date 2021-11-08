export COMPOSE_PROJECT_NAME=monoetna
export NODE_ENV=development

.PHONY: help
help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

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

.PHONY: restart
restart: ## Restarts all projects' processes
	@ make -C docker down
	@ make -C docker up

.PHONY: ps
ps: ## Shows ps of all projects' containers
	@ make -C docker ps

.PHONY: logs
logs: ## Shows logs of all running projects' containers
	@ make -C docker logs

.PHONY: logs-recent
logs-recent: ## For CI
	@ make -C docker logs-recent

.PHONY: update
update: ## Runs all projects' updates, running bundle and npm install
	@ make -C docker update

.PHONY: release
release: ## Builds static docker images staged for release, runs tests against them, and pushes them to dockerhub (requires PUSH_IMAGES=1)
	@ NO_PRUNE=1 $(MAKE) -f Release.mk
	@ docker image prune -f

.PHONE: clean
clean:  ## Cleans many dangling docker references, recovering much disk space.
	docker container prune
	docker images | grep monoetna | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep docker_ | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep metis_ | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep magma_ | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep janus_ | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep timur_ | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep vulcan_ | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker images | grep etnaagent | cut -d ' ' -f1 | xargs -n1 docker image rm || true
	docker image prune
