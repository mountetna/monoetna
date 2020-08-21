export COMPOSE_PROJECT_NAME:=monoetna
COMPOSE_MIXINS:=docker-compose.shared.yml $(COMPOSE_MIXINS)
app_service_name:=${app_name}_app

docker-compose.yml:: $(wildcard ../docker/*.shared.yml) ../docker/default_compose
	COMPOSE_MIXINS="$(COMPOSE_MIXINS)" ../docker/default_compose docker-compose.yml

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
	@ docker-compose down

ps::
	@ docker-compose ps

restart:: docker-ready
	@ docker-compose restart

bash:: docker-ready
	@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm $(app_service_name) bash

logs::
	@ docker-compose logs -f

test:: docker-ready
	@ true

.dockerignore:
	cp ../docker/.dockerignore.template ./.dockerignore

release-build:: .dockerignore
	mkdir -p /tmp/releases
	touch /tmp/releases/success
	rm /tmp/digest
	../docker/build_image -d Dockerfile $(BUILD_REQS) > /tmp/digest
	cat /tmp/digest
	if ! grep "$$(cat /tmp/digest)" /tmp/releases/success && [ -n "$$PULL_IMAGES" ]; then docker pull $(fullTag) || true; fi
	if ! grep "$$(cat /tmp/digest)" /tmp/releases/success; then ../docker/build_image Dockerfile $(BUILD_REQS) -- $(BUILD_ARGS); fi

release::
	make release-build
	if ! grep "$$(cat /tmp/digest)" /tmp/releases/success && ! [ -n "$$NO_TEST" ]; then make release-test; fi
	if  grep "$$(cat /tmp/digest)" /tmp/releases/success && [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi
	cat /tmp/digest >> /tmp/releases/success

release-test::
	true