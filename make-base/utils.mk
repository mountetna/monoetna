.PHONY: always

# Function that find project directories by a file that must exist inside of it
define find_projects_containing
$(sort \
	$(call map,readlink, \
		$(filter-out $(shell dirname $(firstword $(wildcard $(addsuffix /monoetna,.. ../.. ../../..))))/%, \
			$(wildcard \
				$(addsuffix /*/$(1),. .. ../.. ../../..) \
				$(addsuffix /docker/*/$(1),. .. ../.. ../../..) \
				$(addsuffix /swarm/*/$(1),. .. ../.. ../../..) \
			) \
		) \
	) \
)
endef

define find_project_containing
$(shell dirname $(filter %/$(1)/$(2),$(call find_projects_containing,$(2))))
endef

define find_project_file
$(call find_project_containing,$(1),$(2))/$(2)
endef

map = $(foreach a,$(2),$(call $(1),$(a)))
readlink = $(shell readlink -m $(1))
dirname = $(shell basename $(shell dirname $(1)))

define find_updated_targets
$(shell $(call find_project_file,make-base,find-updated-sources) $(1))
endef

define updated_source_files
$(call map,find_updated_targets,$(source_dirs))
endef

define updated_build_files
$(call find_updated_targets,$(call find_project_containing,docker,build_image))
endef

define image_target
$(foreach image,$(shell cat $(1) | echo 'FROM   a' | xargs -n2 echo | cut -d ' ' -f2),include $(call find_project_file,$(image),build.mk))

image.marker: $(updated_source_files) $(updated_build_files) $(call dependent_image_markers,$(1))
	$(call find_project_file,docker,build_image) $(1) $(source_dirs)
	touch $@
endef

define dependent_image_markers
$(foreach image,$(shell cat $(1) | echo 'FROM   a' | xargs -n2 echo | cut -d ' ' -f2),$(call find_project_containing,$(image),build.mk)/image.marker)
endef
