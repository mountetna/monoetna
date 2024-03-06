class Vulcan
  class SSH
    def initialize(net_ssh_instance)
      @ssh = net_ssh_instance
    end
    def invoke_ssh_command(command)
      stdout_data = ""
      stderr_data = ""
      exit_status = nil

      @ssh.open_channel do |channel|
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

    def mkdir(dir)
      # Make project directory if it doesnt exist
      command = "mkdir -p #{dir}"
      invoke_ssh_command(command)
    end

    def dir_exists?(dir)
      command = "[ -d #{dir} ] && echo 'Directory exists.' || echo 'Directory does not exist.'"
      out = invoke_ssh_command(command)
      out[:stdout].chomp == 'Directory exists.'
    end

    def clone(branch, repo, target_dir)
      command = "git clone -b #{branch} #{repo} #{target_dir}"
      invoke_ssh_command(command)
    end


  end
end