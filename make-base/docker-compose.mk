export COMPOSE_PROJECT_NAME:=monoetna
COMPOSE_MIXINS:=docker-compose.shared.yml $(COMPOSE_MIXINS)
app_service_name:=${app_name}_app
baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)
baseFeTag:=${app_name}_app_fe
fullFeTag:=$(IMAGES_PREFIX)$(baseFeTag)$(IMAGES_POSTFIX)
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
	# Build this project's image using the BUILD_REQS and BUILD_ARGS if there is a Dockerfile to build
	if [ -e Dockerfile ]; then ../docker/build_image Dockerfile $(BUILD_REQS) -- $(BUILD_ARGS); fi
	# Build this project's app_fe image iff there is a Dockerfile in the expected place.
	if [ -e $(baseFeTag)/Dockerfile ]; then ../docker/build_image $(baseFeTag)/Dockerfile -- $(BUILD_ARGS); fi

# Create two markers -- one tracks when release tests are run, the other when builds are complete.
# Important is that both of these files are initialized to the same timestamp, so that tests would not
# be invoked unless a successful build_image updates the build marker.
/tmp/etna-build-markers/$(baseTag):
	mkdir -p /tmp/etna-build-markers
	touch -d "197001010101" $<

/tmp/etna-build-markers/$(baseTag).release-test: /tmp/etna-build-markers/$(baseTag)
	touch -d "197001010101" $<

# This step is called and extended to force a run of release image tests and updates the marker that they have been
# executed.
release-test::
	touch /tmp/etna-build-markers/$(baseTag).release-test

# This step invokes the release iff the build marker is newer than the test marker. ie,
# the image has been built since the last test.
.PHONY: release-test-if-stale
release-test-if-stale: /tmp/etna-build-markers/$(baseTag).release-test /tmp/etna-build-markers/$(baseTag)
	if [[ /tmp/etna-build-markers/$(baseTag) -nt /tmp/etna-build-markers/$(baseTag).release-test ]]; then make release-test; fi

release:: /tmp/etna-build-markers/$(baseTag) /tmp/etna-build-markers/$(baseTag).release-test
	make release-build
	if ! [ -n "$$NO_TEST" ]; then make release-test-if-stale; fi
	if [ -e Dockerfile ]; then if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi; fi
	if [ -e $(baseFeTag)/Dockerfile ]; then if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullFeTag); fi; fi

update::
	@ docker-compose run --rm -e FULL_BUILD=1 -e UPDATE_STATE=1 ${app_service_name} echo 'Updated'