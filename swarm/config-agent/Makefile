include ../../make-base/stubs.mk

run-image-test::
	docker run --rm -v /var/run/docker.sock:/var/run/docker.sock:ro -v $$PWD/test:/script/ config-agent