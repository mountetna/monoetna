include ../make-base/stubs.mk
include ../make-base/docker-compose.mk

test::
	@ docker-compose run --rm etna_app bundle exec rspec
	@ docker-compose run -e RUN_NPM_INSTALL=1 -e SKIP_RUBY_SETUP=1 --rm etna_app npm test

bash::
	@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm etna_app bash
