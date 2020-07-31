export COMPOSE_PROJECT_NAME:=monoetna
export COMPOSE_MIXINS:=docker-compose.shared.yml $(COMPOSE_MIXINS)
app_service_name:=${app_name}_app

docker-compose.yml:: $(wildcard ../docker/*.shared.yml) ../docker/default_compose
	../docker/default_compose docker-compose.yml

.PHONY: images
images: docker-compose.yml
	@ make -C ../docker build

config-ready:: docker-compose.yml
	@ true

docker-ready:: images config-ready
	@ true

up:: docker-ready
	@ docker-compose up -d

down:: docker-compose.yml
	@ docker-compose down --remove-orphans

ps::
	@ docker-compose ps

restart:: docker-ready
	@ docker-compose restart

bash:: docker-ready
	@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm $(app_service_name) bash

install-into-volumes:: docker-ready
	@ docker-compose run -e RUN_NPM_INSTALL=1 --rm $(app_service_name) echo 1

logs::
	@ docker-compose logs -f
