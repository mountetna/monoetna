include ../make-base/stubs.mk

app_name=archimedes
include ../make-base/docker-compose.mk

app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
app_service_name:=${app_name}_app
inputs_dir:=$(shell mktemp -d)
outputs_dir:=$(shell mktemp -d)

release-test:: docker-ready
	docker run --rm -e RELEASE_TEST=1 -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} $(fullTag) /entrypoints/development.sh bundle exec rspec

test:: docker-ready
	@ docker-compose run -e ${app_name_capitalized}_ENV=test -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} --rm ${app_service_name} bundle exec rspec
