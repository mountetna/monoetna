# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Metis = Metis.new()

Metis.add_route('GET', '/', 'ClientController#index')

Metis.add_route('POST', '/upload-start', 'UploadController#start_upload')
Metis.add_route('POST', '/upload-blob', 'UploadController#upload_blob')
Metis.add_route('POST', '/upload-pause', 'UploadController#pause_upload')
Metis.add_route('POST', '/upload-stop', 'UploadController#stop_upload')

Metis.add_route('POST', '/magma-end-point', 'MagmaController#magma_end_point')