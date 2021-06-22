class Polyphemus
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path('../views', __dir__)
    def initialize(request, action=nil)
      super
    end

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Polyphemus.instance.config(:token_name),
        timur_host: Polyphemus.instance.config(:timur)&.dig(:host),
        janus_host: Polyphemus.instance.config(:janus)&.dig(:host),
        metis_host: Polyphemus.instance.config(:metis)&.dig(:host),
        vulcan_host: Polyphemus.instance.config(:vulcan)&.dig(:host)
      }.to_json
    end
  end
end
