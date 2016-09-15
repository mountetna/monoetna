# routes.rb

# This file initiates Metis and sets the routing of the http requests
Metis = Metis.new
Metis.add_route('GET', '/', 'Controller#index')
Metis.add_route('POST', '/magma-end-point', 'Controller#magma_end_point')
Metis.add_route('POST', '/upload', 'Controller#upload')