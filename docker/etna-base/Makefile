include ../../make-base/stubs.mk

release-build::
	make -f build.mk $(shell dirname $$(readlink -m $(firstword $(MAKEFILE_LIST))))/image.marker
