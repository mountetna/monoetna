app_service_name:=${app_name}_app
app_db_name:=${app_name}_db
app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
EXTRA_DOCKER_ARGS:=

config.yml: config.yml.template
	$(call find_project_file,build_support,maybe-move-config)

config-ready:: $(app_name)_app_fe config.yml Dockerfile
	@ true

$(app_name)_app_fe:
	cp -r $(call find_project_file,etna-base,app_fe) $(app_name)_app_fe
	ln -s ../$(app_name)/$(app_name)_app_fe ../docker/$(app_name)_app_fe

psql:: docker-ready
	@ docker-compose run -e PGPASSWORD=password --rm ${app_service_name} psql -h ${app_db_name} -U developer -d ${app_name}_development

Dockerfile:
	cp $(call find_project_file,etna-base,Dockerfile.etna-ruby.default) Dockerfile

run-image-test::
	docker-compose up -d $(app_db_name) || true
	docker run --rm $(EXTRA_DOCKER_ARGS) -e $(app_name_capitalized)_ENV=test \
			-e APP_NAME=$(app_name) -e RELEASE_TEST=1 -e CI_SECRET=$${CI_SECRET} \
			-e IS_CI=$${IS_CI} -e WAIT_FOR_DB=1 -e UPDATE_STATE=1 \
			--network monoetna_default $(fullTag) \
			/entrypoints/development.sh bundle exec rspec

update-ready::
	docker-compose up -d $(app_db_name)
