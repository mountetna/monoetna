class Janus
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path('../views', __dir__)

    def initialize(request, action)
      super
      @token = nil
    end

    private

    def success_json(hash = {})
      success(hash.to_json, 'application/json')
    end

    def h(s)
      ERB::Util.html_escape(s)
    end

    # Quick check that the email is in a somewhat valid format.
    def email_valid?(eml)
      eml =~ /\A[\w+\-.]+@[a-z\d\-]+(\.[a-z\d\-]+)*\.[a-z]+\z/
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
