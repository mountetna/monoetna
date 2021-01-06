include ../make-base/stubs.mk

.projects-mark:
	set +e
	git clone git@github.com:mountetna/redcap-projects.git lib/etls/redcap/projects
	@ touch .projects-mark

docker-ready:: .projects-mark
	@ true

app_name=polyphemus
include ../make-base/etna-ruby.mk
include ../make-base/docker-compose.mk
