# client_controller.rb
# This controller serves the client pages and code.

class ClientController < Metis::Controller
  def index
    erb_view :index
  end

  def config_json
    {
      project_name: @params[:project_name],
      token_name: Metis.instance.config(:token_name),
      timur_host: Metis.instance.config(:timur)&.dig(:host),
      vulcan_host: Metis.instance.config(:vulcan)&.dig(:host),
      janus_host: Metis.instance.config(:janus)&.dig(:host),
      polyphemus_host: Metis.instance.config(:polyphemus)&.dig(:host)
    }.to_json
  end
end
