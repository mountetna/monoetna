# client_controller.rb
# This controller serves the client pages and code.

class ClientController < Metis::Controller
  def index
    @token_name = Metis.instance.config(:token_name)
    @janus_host = Metis.instance.config(:auth_redirect)
    erb_view :index
  end
end
