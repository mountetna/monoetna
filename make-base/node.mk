app_service_name=${app_name}_app
test::
				@ docker-compose run -e RUN_NPM_INSTALL=1 -e SKIP_RUBY_SETUP=1 --rm ${app_service_name} npm test

setup-links::
	@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm ${app_service_name} [ -e ../etna ] && npm link ../etna/packages/etna-js
