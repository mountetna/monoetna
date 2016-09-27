# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Metis = Metis.new()

Metis.add_route('GET', '/', 'ClientController#index')

Metis.add_route('POST', '/upload-init', 'UploadController#initialize_upload')
Metis.add_route('POST', '/upload-blob', 'UploadController#upload_blob')
Metis.add_route('POST', '/upload-pause', 'UploadController#pause_upload')
Metis.add_route('POST', '/upload-stop', 'UploadController#stop_upload')
Metis.add_route('POST', '/upload-query', 'UploadController#query_upload')

Metis.add_route('POST', '/magma-end-point', 'UploadController#magma_end_point')