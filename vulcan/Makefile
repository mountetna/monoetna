include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
include ../make-base/etna-ruby.mk
include ../make-base/node.mk
EXTRA_DOCKER_ARGS:=-v /var/run/docker.sock:/var/run/docker.sock:ro -v $${PWD}/../vulcan/spec/data:/app/spec/data -e HOST_DIR_MAP=/app=$${PWD}/../vulcan $(EXTRA_DOCKER_ARGS)

# Ensures that the etnaagent/archimedes image is ready
compose-ready::
	@ IMAGES_PREFIX=etnaagent/ IMAGES_POSTFIX=:latest make -C ../archimedes release-build

start-storybook: docker-ready
	docker run --rm \
	 -e VULCAN_ENV=development \
	 -e APP_NAME=vulcan \
	 -v $$PWD:/app \
	 -v $$PWD/../etna:/etna \
	 -v vulcan-node-modules:/app/node_modules \
	 --entrypoint /entrypoints/development.sh \
	 -p 9700:6006 \
	 -p 6006:6006 \
	 --name vulcan-storybook \
	 etna-base-dev npm run storybook || true

stop-storybook:
	docker stop vulcan-storybook || true