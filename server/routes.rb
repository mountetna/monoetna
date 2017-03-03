# routes.rb
# This file initiates Metis and sets the routing of the http requests.

Polyphemus = Polyphemus.new()

Polyphemus.add_route('GET', '/', 'ClientController#index')