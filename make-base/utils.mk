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
$(call map,dirname,$(filter %/$(1)/$(2),$(call find_projects_containing,$(2))))
endef

define find_project_file
$(addsuffix /$(2),$(call find_project_containing,$(1),$(2)))
endef

map = $(foreach a,$(2),$(call $(1),$(a)))
readlink = $(shell readlink -m $(1))
dirname = $(shell dirname $(1))

define find_updated_sources
$(shell $(call find_project_file,make-base,find-updated-sources) $(1))
endef
define updated_build_files
$(call find_updated_sources,$(call find_project_containing,docker,build_image))
endef

images=$(shell cat $(1) | grep FROM | xargs -n2 echo | cut -d ' ' -f2)

define image_target
$(info $(1) <- $(call images,$(1)))
$(shell if ! test -e $(shell dirname $(1))/image.marker; then touch -t 0001011000 $(shell dirname $(1))/image.marker; fi)
$(addprefix include ,$(foreach image,$(call images,$(1)),$(call find_project_file,$(image),build.mk)))
$(info sources: $(call find_updated_sources,$(shell dirname $(1))))
$(info goal: $(shell dirname $(1)))
$(info deps: $(call find_updated_sources,$(shell dirname $(1))) $(updated_build_files) $(call dependent_image_markers,$(call images,$(1))) $(1))
$(shell dirname $(1))/image.marker: $(call find_updated_sources,$(shell dirname $(1))) $(updated_build_files) $(call dependent_image_markers,$(call images,$(1))) $(1)
	$(call find_project_file,docker,build_image) $(1)
	touch $@
endef

define dependent_image_markers
$(addsuffix /image.marker,$(foreach image,$(1),$(call find_project_containing,$(image),build.mk)))
endef
