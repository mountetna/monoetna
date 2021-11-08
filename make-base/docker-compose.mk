export COMPOSE_PROJECT_NAME:=monoetna
app_service_name:=${app_name}_app
containerSh:=bash

.PHONY: images

ifneq ($(wildcard ./docker-compose.yml),)
$(foreach df,$(call buildable_dependent_compose_image_dockerfiles,$(shell dirname $(firstword $(MAKEFILE_LIST)))/docker-compose.yml),$(eval $(call image_target2,$(shell dirname $(df)),$(call buildable_dependent_image_dockerfiles,$(df)))))
images: config-ready $(foreach df,$(call buildable_dependent_compose_image_dockerfiles,$(shell dirname $(firstword $(MAKEFILE_LIST)))/docker-compose.yml),$(shell dirname $(df))/image.marker)
endif

# docker-compose.yml is a file that must be hand written by the developer in any
# project that uses docker-compose to manage its dev environment.
config-ready:: docker-compose.yml .dockerignore
	@ true

compose-ready:: docker-compose.yml images
	@ true

up:: compose-ready
	@ docker-compose up -d

down:: docker-compose.yml
	@ docker-compose down

ps:: compose-ready
	@ docker-compose ps

restart:: compose-ready
	@ docker-compose restart

bash:: compose-ready
	@ docker-compose run --rm $(app_service_name) $(containerSh)

logs:: compose-ready
	@ docker-compose logs -f

run-image-test:: compose-ready
	@ docker-compose up

.dockerignore:
	cp $(call find_project_file,docker,.dockerignore.template) .dockerignore

update:: compose-ready
	@ docker-compose run --rm -e FULL_BUILD=1 -e UPDATE_STATE=1 ${app_service_name} echo 'Updated'