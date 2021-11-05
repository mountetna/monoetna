include $(shell dirname $(lastword $(MAKEFILE_LIST)))/../../make-base/utils.mk

$(info $(call find_project_file,etna-base,Dockerfile))
$(eval $(call image_target,$(call find_project_file,etna-base,Dockerfile)))
