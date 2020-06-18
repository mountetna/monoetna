# SHELL := /bin/bash
DB_PORT = $(shell docker inspect --format='{{(index (index .NetworkSettings.Ports "5432/tcp") 0).HostPort}}' janus_db_1)

help: ## Display help text
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) /dev/null | \
		sed 's/^[^:]*://' | sort | \
		awk -F':.*?## ' '{printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: help
.DEFAULT_GOAL := help

vendor/bundle: Gemfile Gemfile.lock docker/app/Dockerfile
				@ $(MAKE) bundle
				@ touch vendor/bundle

.docker-build-mark: $(wildcard docker/**/*) docker-compose.yml
				docker-compose rm -f janus_db
				docker-compose pull janus_db
				docker-compose build
				@ touch .docker-build-mark

config.yml:
				cp config.yml.template config.yml

.PHONY: up
up: config.yml .docker-build-mark ## Starts up the database, worker, and webservers of janus in the background.
				@ docker-compose up -d

.PHONY: down
down: ## Ends background janus processes
				@ docker-compose down --remove-orphans

.PHONY: ps
ps: ## Lists status of running janus processes
				@ docker-compose ps

.PHONY: bundle
bundle: ## Executes a bundle install inside of the janus app context.
				docker-compose run --rm janus_app bundle install

.PHONY: build
build: ## Rebuilds the janus docker environment.  Does not clear volumes or databases, just rebuilds code components.
				@ docker-compose build

.PHONY: console
console: ## Starts an irb console inside of the janus app context.
				docker exec -ti "$$(docker ps --format '{{.Names}}' | grep janus_app)" bundle exec irb

.PHONY: migrate
migrate: ## Executes dev and test migrations inside of the janus app context.
				@ docker-compose run --rm janus_app ./bin/janus migrate
				@ docker-compose run -e JANUS_ENV=test --rm janus_app ./bin/janus migrate

.PHONY: test
test: ## Execute (all) rspec tests inside of the janus app context.
				@ docker-compose run -e JANUS_ENV=test --rm janus_app bundle exec rspec

.PHONY: bash
bash: ## Start a bash shell inside of the app context.
				@docker exec -ti "$$(docker ps --format '{{.Names}}' | grep janus_app)" bash

.PHONY: db-port
db-port: ## Print the db port associated with the app.
				@ echo $(DB_PORT)

.PHONY: psql
psql: ## Start a psql shell conntected to the janus development db
				@ docker exec -ti -e PGPASSWORD=password "$$(docker ps --format '{{.Names}}' | grep janus_app)" psql -h janus_db -U developer -d janus_development

.PHONY: logs
logs: ## Follow logs of running containers
				docker-compose logs -f
