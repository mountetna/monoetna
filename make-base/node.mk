app_service_name=${app_name}_app

docker-ready::
	true

run-image-test::
	docker run --rm -e APP_NAME=$(app_name) $(fullTag) npm test
	docker volume rm sync-assets || true
	docker run --rm -v sync-assets:/sync-assets $(fullTag) deploy.sh syncAssets true
	docker volume rm sync-assets || true
