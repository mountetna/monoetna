app_name=deployer
BUILD_REQS=../docker

include ../make-base/stubs.mk
include ../make-base/docker-compose.mk

test::
	docker-compose run deployer_app /test.sh

release-test::
	docker run --rm --mount type=bind,source=/var/run/docker.sock,target=/var/run/docker.sock,readonly --mount type=bind,source=${PWD}/test.sh,target=/test.sh,readonly --mount type=bind,source=${PWD}/TestDockerFile,target=/TestDockerFile,readonly deployer /test.sh
