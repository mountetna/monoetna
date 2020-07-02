include ../make-base/stubs.mk

app_name=metis
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk
include ../make-base/node.mk

setup-links::
	@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm ${app_service_name} npm link ../etna/packages/etna-js
