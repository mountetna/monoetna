# client_controller.rb
# This controller serves the client pages and code. 

class ClientController < Metis::Controller
  def index
    @token_name = Metis.instance.config(:token_name)
    erb_view :index
  end
end
