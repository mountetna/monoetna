baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)
dockerOptions:=-u $(shell id -u):0 -v $$PWD/test_dags:/opt/airflow/dags  -v $$PWD/opt/plugins:/opt/airflow/plugins -v $$PWD/opt/swarm_executor:/opt/airflow/swarm_executor -v /var/run/docker.sock:/var/run/docker.sock:ro

.PHONY: image
image:
	../../docker/build_image -f Dockerfile ../../docker/etna-base-dev
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi

.PHONY: bash
bash: image
	docker run --rm -it $(dockerOptions) $(baseTag) bash

.PHONY: test
test: image
	docker run --rm $(dockerOptions) $(baseTag) python -m pytest

