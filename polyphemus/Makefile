include ../make-base/stubs.mk

.projects-mark:
	git clone git@github.com:mountetna/redcap-projects1.git lib/etls/redcap/projects || clone_failed=1
	if [ ${clone_failed:-0} -eq 1 ]; then \
		echo "Could not clone REDCap projects -- do you have the right permissions?"; \
	fi
	@ touch .projects-mark

docker-ready:: .projects-mark
	@ true

app_name=polyphemus
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk
