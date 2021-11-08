include ../make-base/stubs.mk
include ../make-base/docker-compose.mk
include ../make-base/etna-ruby.mk
include ../make-base/node.mk

release-test::
	# Ensure that the janus css is created
	docker run --rm $(fullTag) ls public/css/janus.bundle.css
