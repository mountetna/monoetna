baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)
dockerOptions:=-u $(shell id -u):0 -v $$PWD/opt/plugins:/opt/airflow/plugins -v /var/run/docker.sock:/var/run/docker.sock:ro

.PHONY: image
image:
	../../docker/build_image -f Dockerfile ../../docker/etna-base-dev
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi

.PHONY: bash
bash: image
	docker run --rm -it $(dockerOptions) $(baseTag) bash

.PHONY: test
test: image
	docker run --rm $(dockerOptions) $(baseTag) -v $$PWD/test_dags:/opt/airflow/dags python -m pytest

