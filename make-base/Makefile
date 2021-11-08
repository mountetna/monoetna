include ./stubs.mk

.PHONY: _test-find-projects
_test-find-projects:
	[[ "18" == "$$(echo $(call find_projects_containing,Makefile) | wc -w)" ]] || exit 1
	[[ "1" == "$$(echo $(filter %etna-base-dev/Makefile,$(call find_projects_containing,Makefile)) | wc -w)" ]] || exit 1
	# Works from root
	[[ "$$(make -f ./stubs.mk _test-find-all)" != "$$(make -C ../ -f make-base/stubs.mk _test-find-all)" ]] || exit 1
	# Works from inside inner directory.
	[[ "$$(make -f ./stubs.mk _test-find-all)" != "$$(make -C ../docker -f ../make-base/stubs.mk _test-find-all)" ]] || exit 1

.PHONY: _test-find-project
_test-find-project:
	[[ "1" == "$$(echo $(filter %etna-base-dev,$(call find_project_containing,etna-base-dev,Makefile)) | wc -w)" ]] || exit 1
	[[ "0" == "$$(echo $(call find_project_containing,etna-base-dev,doesnotexist) | wc -w)" ]] || exit 1
	[[ "0" == "$$(echo $(call find_project_containing,etna-base-d,Makefile) | wc -w)" ]] || exit 1

.PHONY: _test-find-project-file
_test-find-project-file:
	[[ "1" == "$$(echo $(filter %etna-base-dev/Makefile,$(call find_project_file,etna-base-dev,Makefile)) | wc -w)" ]] || exit 1
	[[ "0" == "$$(echo $(call find_project_file,etna-base-dev,doesnotexist) | wc -w)" ]] || exit 1
	[[ "0" == "$$(echo $(call find_project_file,etna-base-d,Makefile) | wc -w)" ]] || exit 1

release-test:: _test-find-projects _test-find-project _test-find-project-file
