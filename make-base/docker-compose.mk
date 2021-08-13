export COMPOSE_PROJECT_NAME:=monoetna
COMPOSE_MIXINS:=docker-compose.shared.yml $(COMPOSE_MIXINS)
app_service_name:=${app_name}_app
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)
containerSh:=bash

docker-compose.yml:: $(wildcard ../docker/*.shared.yml) ../docker/default_compose
	COMPOSE_MIXINS="$(COMPOSE_MIXINS)" ../docker/default_compose docker-compose.yml

.PHONY: images
images: config-ready docker-compose.yml
	@ make -C ../docker build

config-ready:: docker-compose.yml
	@ true

docker-ready:: images
	@ true

up:: docker-ready
	@ docker-compose up -d

down:: docker-compose.yml
	@ docker-compose down

ps::
	@ docker-compose ps

restart:: docker-ready
	@ docker-compose restart

bash:: docker-ready
	@ docker-compose run --rm $(app_service_name) $(containerSh)

logs::
	@ docker-compose logs -f

test:: docker-ready
	@ true

.dockerignore:
	if [ -e Dockerfile]; then cp ../docker/.dockerignore.template ./.dockerignore; fi

release-build:: .dockerignore
	if [ -e Dockerfile ]; then touch /tmp/.release-test; ../docker/build_image Dockerfile $(BUILD_REQS) -- $(BUILD_ARGS); else touch /tmp/etna-build-markers/$(baseTag); fi

/tmp/.release-test: /tmp/etna-build-markers/$(baseTag)
	make release-test
	touch /tmp/.release-test

release::
	mkdir -p /tmp/etna-build-markers
	touch /tmp/etna-build-markers/$(baseTag)
	touch /tmp/.release-test
	make release-build
	if ! [ -n "$$NO_TEST" ]; then make /tmp/.release-test; fi
	if [ -e Dockerfile ]; then if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi; fi

update::
	@ docker-compose run --rm -e FULL_BUILD=1 -e UPDATE_STATE=1 ${app_service_name} echo 'Updated'

release-test::
	true