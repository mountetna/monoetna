export COMPOSE_PROJECT_NAME=monoetna
app_service_name=${app_name}_app

prepare-compose::
	export COMPOSE_MIXINS="docker-compose.shared.yml $COMPOSE_MIXINS"

docker-compose.yml:: prepare-compose $(wildcard ../docker/*.shared.yml) ../docker/default-compose
	../docker/default-compose docker-compose.yml

.PHONY: images
images: docker-compose.yml
	../docker/build_images_of docker-compose.yml

docker-ready:: docker-compose.yml images

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
