app_service_name=${app_name}_app

docker-ready::
	true

release-test::
	docker run --rm $(fullTag) npm test
	docker volume rm sync-assets || true
	docker run --rm -v sync-assets:/sync-assets $(fullTag) deploy syncAssets true
	docker volume rm sync-assets || true
