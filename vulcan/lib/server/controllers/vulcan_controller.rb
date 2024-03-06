require 'etna'

class Vulcan
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def initialize(request, action = nil)
      super
      escaped_params
    end

    def invoke_ssh_command(command)
      stdout_data = ""
      stderr_data = ""
      exit_status = nil

      Vulcan.instance.ssh.open_channel do |channel|
        channel.exec(command) do |ch, success|
          unless success
            raise "Command execution failed: #{command}"
          end

          channel.on_data do |_, data|
            stdout_data += data
          end

          channel.on_extended_data do |_, type, data|
            # Note: here stderr is used for both errors and information messages
            stderr_data += data
          end

          channel.on_request("exit-status") do |_, data|
            exit_status = data.read_long
          end
        end
      end

      Vulcan.instance.ssh.loop

      if exit_status != 0
        raise "Command exited with status #{exit_status}: #{command}"
      end

      {command: command, stdout: stdout_data, stderr_or_info: stderr_data, exit_status: exit_status }
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
