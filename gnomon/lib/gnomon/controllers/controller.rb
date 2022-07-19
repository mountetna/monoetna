class Gnomon
  class Controller < Etna::Controller
    def config_json
      {
        project_name: @params[:project_name],
        token_name: Gnomon.instance.config(:token_name)
      }.merge(config_hosts).to_json
    end
  end
end
