app_name=etna
include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
BUILD_ARGS:=--build-arg APP_NAME=$(app_name) $(BUILD_ARGS)
export BUILD_REQS:=../docker/etna-base $(BUILD_REQS)

bash::
	@ docker-compose run --rm etna_app bash

release-test::
	docker network create monoetna_default
	docker run --rm -e APP_NAME=etna -e RELEASE_TEST=1 -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} --network monoetna_default $(fullTag) /entrypoints/development.sh bundle exec rspec
	docker run --rm -e APP_NAME=etna -e RELEASE_TEST=1 -e IS_CI=$${IS_CI} --network monoetna_default $(fullTag) /entrypoints/development.sh npm test