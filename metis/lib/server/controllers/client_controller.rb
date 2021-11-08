# client_controller.rb
# This controller serves the client pages and code.

class ClientController < Metis::Controller
  def index
    erb_view :index
  end

  def config_json
    {
      project_name: @params[:project_name],
      token_name: Metis.instance.config(:token_name)
    }.merge(config_hosts).to_json
  end
end
