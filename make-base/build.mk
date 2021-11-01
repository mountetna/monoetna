.PHONY: always

define make_build_targets
$(1)/build/image: $(1)/build/latest-change
	mkdir -p $(1)/build

$(1)/build/latest-change: always
	mkdir -p $(1)/build

.PHONY: $(1)/clean
$(1)/clean:
	rm -rf $(1)/build
endef

$(eval $(call make_build_targets,\
		$(shell dirname $(abspath $(lastword $(MAKEFILE_LIST)))),\
		$(shell dirname $(shell readlink -m $(lastword $(MAKEFILE_LIST))))\
		)\
)