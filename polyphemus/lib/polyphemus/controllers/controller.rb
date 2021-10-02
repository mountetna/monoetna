class Polyphemus
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path('../views', __dir__)
    def initialize(request, action=nil)
      super
    end

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Polyphemus.instance.config(:token_name)
      }.merge(config_hosts).to_json
    end
  end
end
