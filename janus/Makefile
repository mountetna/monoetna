BUILD_ARGS:=--build-arg RUN_NPM_INSTALL=1 $(BUILD_ARGS)

include ../make-base/stubs.mk

app_name=janus
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk

release-test::
	# Ensure that the janus css is created
	docker run --rm $(fullTag) ls lib/client/css/janus.css