# client_controller.rb
# This controller serves the client pages and code. 

class ClientController < Metis::Controller
  def index
    view :index
  end

  def user
    view :user
  end
end
