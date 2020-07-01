export COMPOSE_PROJECT_NAME=monoetna
app_service_name=${app_name}_app

.docker-build-mark: $(wildcard docker/**/*) docker-compose.yml
				@ $(MAKE) build

up:: .docker-build-mark
				@ docker-compose up -d

down::
				@ docker-compose down --remove-orphans

ps::
				@ docker-compose ps

build::
				@ docker-compose build
				@ touch .docker-build-mark

restart::
				@ docker-compose restart

bash::
				@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm $(app_service_name) bash

logs::
				@ docker-compose logs -f
