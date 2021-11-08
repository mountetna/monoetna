include make-base/utils.mk

define setup_all_releases
$(foreach df,$(call find_projects_containing,Dockerfile),$(call image_target2,$(shell dirname $(df)),$(call buildable_dependent_image_dockerfiles,$(df))))
$(foreach df,$(call find_projects_containing,Dockerfile),$(call release_image_target1,$(shell dirname $(df))))
$(foreach df,$(call find_projects_containing,Dockerfile),$(call test_image_target1,$(shell dirname $(df))))
endef

release: $(foreach df,$(call find_projects_containing,Dockerfile),$(shell dirname $(df))/image-release.marker)
$(info release: $(foreach df,$(call find_projects_containing,Dockerfile),$(shell dirname $(df))/image-release.marker))

$(eval $(call setup_all_releases))
