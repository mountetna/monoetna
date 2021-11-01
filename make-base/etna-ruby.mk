app_service_name:=${app_name}_app
app_db_name:=${app_name}_db
app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
BUILD_ARGS:=--build-arg APP_NAME=$(app_name) $(BUILD_ARGS)
EXTRA_DOCKER_ARGS:=
export BUILD_REQS:=../docker/etna-base $(BUILD_REQS)

config.yml: config.yml.template
	../make-base/maybe-move-config

config-ready:: config.yml
	@ true

irb::
	@ docker-compose run --rm $app_service_name bundle exec irb

test:: docker-ready
	@ docker-compose run -e ${app_name_capitalized}_ENV=test -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} --rm ${app_service_name} bundle exec rspec

psql:: docker-ready
	@ docker-compose run -e -e PGPASSWORD=password --rm ${app_service_name} psql -h ${app_db_name} -U developer -d ${app_name}_development

Dockerfile:
	cp ../docker/etna-base/release/Dockerfile .

release-build:: .dockerignore

release-test:: docker-ready
	docker-compose up -d $(app_db_name)
	docker run --rm $(EXTRA_DOCKER_ARGS) -e $(app_name_capitalized)_ENV=test -e APP_NAME=$(app_name) -e RELEASE_TEST=1 -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} -e WAIT_FOR_DB=1 -e UPDATE_STATE=1 --network monoetna_default $(fullTag) /entrypoints/development.sh bundle exec rspec
