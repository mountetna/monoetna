app_service_name=${app_name}_app

docker-ready::
	true

release-test::
	docker run --rm $(fullTag) npm test
