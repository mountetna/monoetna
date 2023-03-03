include ../make-base/stubs.mk
include ../make-base/docker-compose.mk

app_name_capitalized:=$(shell echo ${app_name} | tr [a-z] [A-Z])
app_service_name:=${app_name}_app
fullTag:=$(IMAGES_PREFIX)$(app_name)$(IMAGES_POSTFIX)
