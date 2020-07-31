app_service_name=${app_name}_app
app_db_name=${app_name}_db
app_name_capitalized=$(shell echo ${app_name} | tr [a-z] [A-Z])

prepare-compose::
	export COMPOSE_MIXINS="docker-compose.etna-app.shared.yml $COMPOSE_MIXINS"

config.yml: config.yml.template
	@ [[ ! -e config.yml || (read -p "$DIR/config.yml.template has been updated recently, copy it into config.yml? y/n" && $REPLY =~ [yY]) ]] && cp config.yml.template config.yml

config-ready:: config.yml
	@ true

irb::
	@ docker-compose run --rm $app_service_name bundle exec irb

migrate::
	@ docker-compose run --rm ${app_service_name} ./bin/${app_name} migrate
	@ docker-compose run -e ${app_name_capitalized}_ENV=test --rm ${app_service_name} ./bin/${app_name} migrate

test:: docker-ready
	@ docker-compose run -e ${app_name_capitalized}_ENV=test --rm ${app_service_name} bundle exec rspec

psql:: docker-ready
	@ docker-compose run -e SKIP_RUBY_SETUP=1 -e PGPASSWORD=password --rm ${app_service_name} psql -h ${app_db_name} -U developer -d ${app_name}_development
