include ../make-base/stubs.mk

.projects-mark:
	git clone git@github.com:mountetna/redcap-projects.git lib/etls/redcap/projects || REDCAP_CLONE_FAILED=1
	git clone git@github.com:mountetna/renaming-projects.git lib/etls/renaming/projects || RENAMING_CLONE_FAILED=1
	if [ -n "$(REDCAP_CLONE_FAILED)" ]; then \
		echo "Could not clone REDCap projects -- do you have the right permissions?"; \
	elif [ -n "$(RENAMING_CLONE_FAILED)" ] ; then \
		echo "Could not clone Renaming projects -- do you have the right permissions?"; \
	else \
		touch .projects-mark; \
	fi
	@ true

docker-ready:: .projects-mark
	@ true

app_name=polyphemus
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk
include ../make-base/node.mk
