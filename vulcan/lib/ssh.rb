require 'shellwords'

class Vulcan
  class SSH
    def initialize(net_ssh_instance)
      @ssh = net_ssh_instance
    end

    def mkdir(dir)
      # Make project directory if it doesnt exist
      command = Shellwords.join(["mkdir", "-p", dir])
      invoke_ssh_command(command)
    end

    def rmdir(dir, allowed_directories)
      # Check if the dir is in the allowed_directories
      if allowed_directories.any? { |allowed_dir| dir.start_with?(allowed_dir) }
        command = Shellwords.join(["rm", "-r", "-f", dir])
        invoke_ssh_command(command)
      else
        raise ArgumentError, "Directory #{dir} is not in the list of allowed directories."
      end
    end

    def dir_exists?(dir)
      command = "[ -d #{dir} ] && echo 'Directory exists.' || echo 'Directory does not exist.'"
      out = invoke_ssh_command(command)
      out[:stdout].chomp == 'Directory exists.'
    end

    def clone(repo, branch, target_dir)
      command = Shellwords.join(['git', 'clone', '-b', branch, repo, target_dir])
      invoke_ssh_command(command)
    end

    private
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

  end
end