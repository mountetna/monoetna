require 'etna'
require 'shellwords'

class Vulcan
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def initialize(request, action = nil)
      super
      escaped_params
      transform_array_values
    end

    private

    def task_token
      janus_client = Etna::Clients::Janus.new(token: @user.token, host: Vulcan.instance.config(:janus)[:host])
      janus_client.generate_token("task", project_name: @params[:project_name], read_only: true)
    end

    def escaped_params
      @escaped_params ||= @params.transform_values { |value| Shellwords.escape(value.to_s) }
    end

    def transform_array_values
      @params.each do |key, value|
        if value.is_a?(Array) && value.all?(&:empty?)
          @params[key] = []
        end
      end
    end

    def config_json
      {
        project_name: @params[:project_name],
        token_name: Vulcan.instance.config(:token_name)
      }.merge(config_hosts).to_json
    end

    def token
      @token ||= @request.cookies[Vulcan.instance.config(:token_name)]
    end

    # This recursively converts string booleans ("true"/"false") into actual boolean values.
    def to_valid_json(obj_to_transform)
      JSON.generate(recursive_transform(obj_to_transform))
    end

    # Helper method to recursively transform hash values.
    def recursive_transform(value)
      case value
      when Hash
        value.transform_values { |v| recursive_transform(v) }
      when Array
        value.map { |v| recursive_transform(v) }
      when String
        normalized = value.strip.downcase
        return true if normalized == 'true'
        return false if normalized == 'false'
        value
      else
        value
      end
    end
  end
end
