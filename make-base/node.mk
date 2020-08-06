app_service_name=${app_name}_app
BUILD_ARGS:=--build-arg RUN_NPM_INSTALL=1 $(BUILD_ARGS)
BUILD_REQS:=../etna $(BUILD_REQS)

docker-ready::
	true

test:: docker-ready
	@ docker-compose run -e RUN_NPM_INSTALL=1 -e SKIP_RUBY_SETUP=1 --rm ${app_service_name} npm test

release-test::
	docker run --rm $(fullTag) npm test
