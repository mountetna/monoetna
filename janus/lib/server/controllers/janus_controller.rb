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
      }.merge(config_hosts).to_json
    end
  end
end
