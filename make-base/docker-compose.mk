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
	@ docker-compose run -e SKIP_RUBY_SETUP=1 -e SKIP_PYTHON_SETUP=1 --rm $(app_service_name) $(containerSh)

logs::
	@ docker-compose logs -f

test:: docker-ready
	@ true

.dockerignore:
	test -e Dockerfile && cp ../docker/.dockerignore.template ./.dockerignore

release-build:: .dockerignore
	test -e Dockerfile && ../docker/build_image Dockerfile $(BUILD_REQS) -- $(BUILD_ARGS)

release::
	make release-build
	if ! [ -n "$$NO_TEST" ]; then make release-test; fi
	test -e Dockerfile && if   [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi

release-test::
	true