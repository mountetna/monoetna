include ../make-base/utils.mk

$(eval $(call image_target,$(call find_project_file,etna,Dockerfile)))
