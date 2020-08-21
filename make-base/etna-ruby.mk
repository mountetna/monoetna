app_service_name:=${app_name}_app
app_db_name:=${app_name}_db
app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
COMPOSE_MIXINS:=docker-compose.etna-app.shared.yml $(COMPOSE_MIXINS)
BUILD_ARGS:=--build-arg SKIP_RUBY_SETUP= $(BUILD_ARGS)

config.yml: config.yml.template
	../make-base/maybe-move-config

config-ready:: config.yml
	@ true

irb::
	@ docker-compose run --rm $app_service_name bundle exec irb

migrate::
	@ docker-compose run --rm ${app_service_name} ./bin/${app_name} migrate
	@ docker-compose run -e ${app_name_capitalized}_ENV=test --rm ${app_service_name} ./bin/${app_name} migrate

test:: docker-ready
	@ docker-compose run -e ${app_name_capitalized}_ENV=test --rm ${app_service_name} bundle exec rspec

psql:: docker-ready
	@ docker-compose run -e SKIP_RUBY_SETUP=1 -e PGPASSWORD=password --rm ${app_service_name} psql -h ${app_db_name} -U developer -d ${app_name}_development

Dockerfile:
	cp ../docker/etna-base/release/Dockerfile .

release:: Dockerfile
	@ true

release-test:: docker-ready
	docker-compose up -d $(app_db_name)
	docker run --rm -e $(app_name_capitalized)_ENV=test -e APP_NAME=$(app_name) -e RELEASE_TEST=1 --network monoetna_default $(fullTag) /entrypoints/development.sh bundle exec rspec
