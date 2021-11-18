include ../../make-base/stubs.mk

baseTag:=$(shell basename "$$(pwd)")
dockerOptions:=-u $(shell id -u):0 -v $$PWD/opt/plugins:/opt/airflow/plugins -v /var/run/docker.sock:/var/run/docker.sock:ro

bash:: release-build
	docker run --rm -it $(dockerOptions) $(baseTag) bash

run-image-test:: release-build
	# docker run --rm $(dockerOptions) -v $$PWD/test_dags:/opt/airflow/dags $(baseTag) python -m pytest
	@ true

