baseTag:=$(shell basename "$$(pwd)")
fullTag:=$(IMAGES_PREFIX)$(baseTag)$(IMAGES_POSTFIX)

.PHONY: image
image:
	../../docker/build_image -f Dockerfile ../../docker/etna-base-dev
	if [ -n "$$PUSH_IMAGES" ]; then docker push $(fullTag); fi

.PHONY: bash
bash: image
	docker run --rm -it -u 0:0 -v $$PWD/opt/plugins:/opt/airflow/plugins -v $$PWD/opt/swarm_executor:/opt/airflow/swarm_executor -v /var/run/docker.sock:/var/run/docker.sock:ro $(baseTag) bash

.PHONY: test
test: image
	docker run --rm -it -u 0:0 -v $$PWD/opt/plugins:/opt/airflow/plugins -v $$PWD/opt/swarm_executor:/opt/airflow/swarm_executor -v /var/run/docker.sock:/var/run/docker.sock:ro $(baseTag) python -m pytest

