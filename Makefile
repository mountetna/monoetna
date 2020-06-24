# SHELL := /bin/bash
export COMPOSE_PROJECT_NAME=monoetna

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

.docker-build-mark: $(wildcard docker/**/*) docker-compose.yml
				docker-compose build
				@ touch .docker-build-mark

.PHONY: up
up: .docker-build-mark ## Stub: etna currently doesn't really run anything itself.
				@ echo etna has no running processes

.PHONY: down
down: ## Stub: etna currently doesn't really run anything itself.
				@ echo etna has no running processes

.PHONY: ps
ps: ## Lists status of running processes
				@ echo etna has no running processes

.PHONY: bundle
bundle: ## Executes a bundle install inside of the etna app context.
				docker-compose run --rm etna_app bundle install

.PHONY: build
build: ## Rebuilds the etna docker environment.  Does not clear npm or bundle caches, just rebuilds code components.
				@ docker-compose build

.PHONY: irb
irb: ## Starts an irb console inside of the running etna app context.
				@ docker-compose run --rm etna_app bundle exec irb

.PHONY: node
node: ## Starts an node console inside of the running etna app context.
				@ docker-compose run --rm etna_app node

.PHONY: test
test: ## Execute (all) rspec and npm tests inside of the etna app context.
				@ docker-compose run --rm etna_app bundle exec rspec
				@ docker-compose run -e RUN_NPM_INSTALL=1 -e SKIP_RUBY_SETUP=1 --rm etna_app npm test

.PHONY: bash
bash: ## Start a bash shell inside of the app context.
				@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm etna_app bash
