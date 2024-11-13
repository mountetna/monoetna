require 'etna'
require 'shellwords'

class Vulcan
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def initialize(request, action = nil)
      super
      escaped_params
    end

    private

    def escaped_params
      @escaped_params ||= @params.transform_values { |value| Shellwords.escape(value.to_s) }
    end

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Vulcan.instance.config(:token_name)
      }.merge(config_hosts).to_json
    end

    def storage
      @storage ||= Vulcan::Storage.new
    end

    def token
      @token ||= @request.cookies[Vulcan.instance.config(:token_name)]
    end
  end
end
