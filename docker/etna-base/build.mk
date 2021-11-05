.PHONY: always
include ../../make-base/utils.mk

etna-base-dev:=$(call find_project_containing,etna-base-dev,build.mk)
include $(etna-base-dev)/build.mk

source_dirs:=$(call find_project_containing,etna,Dockerfile) .

previous.id.marker: always


$(eval $(call image_target,$(etna-base-dev)/image.marker))
image.marker: previous.id.marker $(updated_source_files) $(updated_build_files) $(etna-base-dev)/image.marker
	if [[ "$$($(call find_project_file,docker,build_image) Dockerfile $(source_dirs))" =~ updated ]]; then \
	  touch $@ \
	else
	  if
	fi
