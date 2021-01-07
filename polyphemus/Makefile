include ../make-base/stubs.mk

.projects-mark:
	git clone git@github.com:mountetna/redcap-projects.git lib/etls/redcap/projects || CLONE_FAILED=1
	if [ -z "$(CLONE_FAILED)" ]; then \
		echo "Could not clone REDCap projects -- do you have the right permissions?"; \
	fi
	@ touch .projects-mark

docker-ready:: .projects-mark
	@ true

app_name=polyphemus
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk
