# SHELL := /bin/bash
export COMPOSE_PROJECT_NAME=monoetna

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

.PHONY: restart
restart: ## Restarts background apache process
				@ docker-compose restart

.PHONY: up
up: ## Brings up background apache process
				@ docker-compose up -d

.PHONY: down
down: ## Stops the background apache process
				@ docker-compose down

.PHONY: logs
logs: ## Follow logs of running containers
				docker-compose logs -f
