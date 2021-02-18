class Janus
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path('../views', __dir__)

    private

    def success_json(hash = {})
      success(hash.to_json, 'application/json')
    end

    def h(s)
      ERB::Util.html_escape(s)
    end

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Janus.instance.config(:token_name),
        timur_host: Janus.instance.config(:timur)&.dig(:host),
        metis_host: Janus.instance.config(:metis)&.dig(:host),
      }.to_json
    end
  end
end
