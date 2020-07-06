app_service_name=${app_name}_app
app_db_name=${app_name}_db
app_name_capitalized=$(shell echo ${app_name} | tr [a-z] [A-Z])

config.yml:
				cp config.yml.template config.yml

up:: config.yml
				@ true

irb::
				@ docker-compose run --rm $app_service_name bundle exec irb

migrate::
				@ docker-compose run --rm ${app_service_name} ./bin/${app_name} migrate
				@ docker-compose run -e ${app_name_capitalized}_ENV=test --rm ${app_service_name} ./bin/${app_name} migrate

test:: config.yml
				@ docker-compose run -e ${app_name_capitalized}_ENV=test --rm ${app_service_name} bundle exec rspec

psql::
				@ docker-compose run -e SKIP_RUBY_SETUP=1 --rm ${app_service_name} psql -h ${app_db_name} -U developer -d ${app_name}_development
