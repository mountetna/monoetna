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
