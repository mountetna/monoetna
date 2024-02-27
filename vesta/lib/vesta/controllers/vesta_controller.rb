class Vesta
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path('../views', __dir__)

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Vesta.instance.config(:token_name)
      }.merge(config_hosts).to_json
    end
  end
end
