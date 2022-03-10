include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
app_service_name:=airflow_ui

baseTag:=$(shell basename "$$(pwd)")
pwd:=$(shell echo "$$(pwd)")
# dockerOptions:=-u $(shell id -u):0 -v $$PWD/opt/plugins:/opt/airflow/plugins -v /var/run/docker.sock:/var/run/docker.sock:ro

run-image-test::
	docker run --rm -e USE_MOCKS=1 $(baseTag) bash -c 'pytest'

black:
	docker run -u 0:0 -v $(pwd)/opt/providers:/opt/airflow/providers --rm $(baseTag) bash -c '$$@' -- black /opt/airflow/providers/etna

pdoc: $(shell readlink -m $(shell dirname $(firstword $(MAKEFILE_LIST))))/image-release.marker
	docker run -u 0:0 -v $(pwd)/opt/providers:/opt/airflow/providers --rm $(baseTag) bash -c '$$@' -- pdoc -o /opt/airflow/providers/etna/docs etna
