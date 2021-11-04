include ../../make-base/utils.mk
source_dirs:=$(call find_project_containing,etna,Dockerfile) .


image.target.mk: $(updated_source_files) \
		  $(updated_build_files) \
		  $(call find_project_file,etna-base-dev,image.target.mk)
	$(call find_project_file,docker,build_image) Dockerfile $(source_dirs)
