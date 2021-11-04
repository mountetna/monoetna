export COMPOSE_PROJECT_NAME:=monoetna
app_service_name:=${app_name}_app
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)
baseFeTag:=${app_name}_app_fe
fullFeTag:=$(IMAGES_PREFIX)$(baseFeTag)$(IMAGES_POSTFIX)
containerSh:=bash

.PHONY: images
images: config-ready docker-compose.yml
	@ make -C ../docker build

config-ready:: docker-compose.yml
	@ true

compose-ready:: images
	@ true

up:: compose-ready
	@ docker-compose up -d

down:: docker-compose.yml
	@ docker-compose down

ps::
	@ docker-compose ps

restart:: compose-ready
	@ docker-compose restart

bash:: compose-ready
	@ docker-compose run --rm $(app_service_name) $(containerSh)

logs::
	@ docker-compose logs -f

test:: compose-ready
	@ true

.dockerignore:
	if [ -e Dockerfile]; then cp ../docker/.dockerignore.template ./.dockerignore; fi

release-build:: .dockerignore
	if [ -e Dockerfile ]; then touch /tmp/.release-test; ../docker/build_image Dockerfile $(BUILD_REQS) -- $(BUILD_ARGS); else touch /tmp/etna-build-markers/$(baseTag); fi
	if [ -e $(baseFeTag)/Dockerfile ]; then ../docker/build_image $(baseFeTag)/Dockerfile -- $(BUILD_ARGS); fi

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
	if [ -e $(baseFeTag)/Dockerfile ]; then if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullFeTag); fi; fi

update::
	@ docker-compose run --rm -e FULL_BUILD=1 -e UPDATE_STATE=1 ${app_service_name} echo 'Updated'

release-test::
	true