# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Metis = Metis.new()

Metis.add_route('GET', '/', 'ClientController#index')

Metis.add_route('POST', '/login', 'LogController#log_in')
Metis.add_route('POST', '/logout', 'LogController#log_out')
Metis.add_route('POST', '/check', 'LogController#check_log')

Metis.add_route('POST', '/upload-authorize', 'UploadController#authorize_upload')
Metis.add_route('POST', '/upload-start', 'UploadController#start_upload')
Metis.add_route('POST', '/upload-blob', 'UploadController#upload_blob')
Metis.add_route('POST', '/upload-pause', 'UploadController#pause_upload')
Metis.add_route('POST', '/upload-stop', 'UploadController#stop_upload')

Metis.add_route('POST', '/magma-end-point', 'MagmaController#magma_end_point')