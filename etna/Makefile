include ../make-base/stubs.mk
include ../make-base/docker-compose.mk

bash::
	@ docker-compose run --rm etna_app bash

run-image-test::
	docker run --rm -e FULL_BUILD=1 -e APP_NAME=etna -e CI_SECRET=$${CI_SECRET} -e IS_CI=$${IS_CI} --network monoetna_default $(fullTag) /entrypoints/development.sh bundle exec rspec
	docker run --rm -e FULL_BUILD=1 -e APP_NAME=etna -e IS_CI=$${IS_CI} --network monoetna_default $(fullTag) /entrypoints/development.sh npm test

update::
	@ docker-compose run --rm -e FULL_BUILD=1 -e UPDATE_STATE=1 etna_app bundle install
	@ docker-compose run --rm -e FULL_BUILD=1 -e UPDATE_STATE=1 etna_app bash -c 'cd packages/etna-js && npm install'
