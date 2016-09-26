# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Metis = Metis.new()
Metis.add_route('GET', '/', 'Controller#index')
Metis.add_route('POST', '/magma-end-point', 'Controller#magma_end_point')

Metis.add_route('POST', '/upload-init', 'Controller#initialize_upload')
Metis.add_route('POST', '/upload-blob', 'Controller#upload_blob')
Metis.add_route('POST', '/upload-pause', 'Controller#pause_upload')
Metis.add_route('POST', '/upload-stop', 'Controller#stop_upload')
Metis.add_route('POST', '/upload-query', 'Controller#query_upload')
