require 'etna'

class Vulcan
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    private

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Vulcan.instance.config(:token_name)
      }.merge(config_hosts).to_json
    end

    def storage
      @storage ||= Vulcan::Storage.new
    end

    def redirect_to(path)
      @response.redirect(path,302)
      @response.finish
    end

    def token
      @token ||= @request.cookies[Vulcan.instance.config(:token_name)]
    end
  end
end
