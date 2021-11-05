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
$(shell $(call find_project_file,build_support,find-updated-sources) $(1))
endef
define updated_build_files
$(call find_updated_sources,$(call find_project_containing,build_support,build_image))
endef

images=$(shell cat $(1) | grep FROM | xargs -n2 echo | cut -d ' ' -f2)

# ( dockerfile )
define buildable_dependent_image_dockerfiles
$(foreach image,$(call images,$(1)),$(call find_project_file,$(image),Dockerfile))
endef

# (dockerfile folder, dependent image dockerfiles)
define image_target2

$(shell if ! test -e $(1)/image.marker; then touch -t 0001011000 $(1)/image.marker; fi)
$(foreach df,$(2),$(call image_target2,$(shell dirname $(df)),$(call buildable_dependent_image_dockerfiles,$(df))))

$(1)/image.marker: $(call find_updated_sources,$(1)) $(updated_build_files) $(foreach df,$(2),$(shell dirname $(df))/image.marker)
	$(call find_project_file,build_support,build_image) $(1)/Dockerfile
	touch $$@
endef

define image_target1
$(call image_target2,$(1),$(call buildable_dependent_image_dockerfiles,$(1)/Dockerfile))
endef

define image_target_here
$(call image_target1,$(shell dirname $$(readlink -m $(firstword $(MAKEFILE_LIST)))))
release-build:: $(shell dirname $$(readlink -m $(firstword $(MAKEFILE_LIST))))/image.marker
endef