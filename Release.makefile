include make-base/utils.mk

.PHONY: release
$(foreach df,$(call find_projects_containing,Dockerfile),$(eval $(call image_target2,$(shell dirname $(df)),$(call buildable_dependent_image_dockerfiles,$(df)))))
$(foreach df,$(call find_projects_containing,Dockerfile),$(eval $(call release_image_target1,$(shell dirname $(df)))))
$(foreach df,$(call find_projects_containing,Dockerfile),$(eval $(call test_image_target1,$(shell dirname $(df)))))
release: $(foreach df,$(call find_projects_containing,Dockerfile),$(shell dirname $(df))/image-release.marker)
$(info release: $(foreach df,$(call find_projects_containing,Dockerfile),$(shell dirname $(df))/image-release.marker))
