include ../../make-base/utils.mk
source_dirs:=$(find_project_file development-certs,certs) .


#image.target.mk: $(updated_source_files) \
#		  $(updated_build_files)
#	$(call find_project_file,docker,build_image) Dockerfile $(source_dirs)
