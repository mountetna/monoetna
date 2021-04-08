app_name=etna
include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
BUILD_ARGS:=--build-arg SKIP_RUBY_SETUP= --build-arg APP_NAME=$(app_name) $(BUILD_ARGS)
export BUILD_REQS:=../docker/etna-base $(BUILD_REQS)

test::
	@ docker-compose run --rm etna_app bundle exec rspec
	@ docker-compose run -e RUN_NPM_INSTALL=1 -e SKIP_RUBY_SETUP=1 --rm etna_app npm test
	# We will have to add in a command here to test the R stuff, too.

bash::
	@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm etna_app bash

release-test::
	docker network create monoetna_default
	docker run --rm -e APP_NAME=etna -e RELEASE_TEST=1 -e SKIP_DB_WAIT=1 -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} --network monoetna_default $(fullTag) /entrypoints/development.sh bundle exec rspec
	docker run --rm -e APP_NAME=etna -e RELEASE_TEST=1 -e SKIP_RUBY_SETUP=1 -e RUN_NPM_INSTALL=1 -e SKIP_DB_WAIT=1 --network monoetna_default $(fullTag) /entrypoints/development.sh npm test
	# We will have to add in a command here to test the R stuff, too.