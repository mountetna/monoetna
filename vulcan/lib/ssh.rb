require 'shellwords'

class Vulcan
  class SSH
    def initialize(net_ssh_instance)
      @ssh = net_ssh_instance
    end

    def run_snakemake(dir, command)
      command = "cd #{Shellwords.escape(dir)} && #{command}"
      get_slurm_run_uuid(command, 5)
    end

    def get_snakemake_log(dir, slurm_run_uuid)
      command = "cd #{Shellwords.escape(dir)} && grep -rl #{slurm_run_uuid} .snakemake/log/"
      out = invoke_ssh_command(command)
      out[:stdout].gsub("\n", "")
    end

    def write_file(remote_file_path, content)
      command = "mkdir -p $(dirname #{remote_file_path}) && echo #{Shellwords.escape(content)} > #{remote_file_path}"
      invoke_ssh_command(command)
    end

    def read_yaml_file(remote_file_path)
      command = Shellwords.join(['cat', remote_file_path])
      result = invoke_ssh_command(command)
      YAML.safe_load(result[:stdout])
    end

    def upload_dir(local_dir, remote_dir, recurse)
      @ssh.scp.upload!(local_dir, remote_dir, recursive: recurse)
    end

    def mkdir(dir)
      # Make project directory if it doesnt exist
      command = Shellwords.join(["mkdir", "-p", dir])
      invoke_ssh_command(command)
    end

    def list_dirs(dir)
      command = "ls -d #{Shellwords.escape(dir)}/*/"
      result = invoke_ssh_command(command)
      result[:stdout].split("\n").map { |path| File.basename(path) }
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

    def checkout_tag(dir, tag)
      command = "cd #{Shellwords.escape(dir)} && git checkout tags/#{Shellwords.escape(tag)}"
      invoke_ssh_command(command)
    end

    def get_repo_remote_url(repo_dir)
      command = "cd #{Shellwords.escape(repo_dir)} && git config --get remote.origin.url"
      result = invoke_ssh_command(command)
      result[:stdout].chomp
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
            stderr_data += data
          end

          channel.on_request("exit-status") do |_, data|
            exit_status = data.read_long
          end
        end
      end

      @ssh.loop

      if exit_status != 0
        raise "Command exited with status #{exit_status}. \n Command: #{command} \n Msg: #{stderr_data}"
      end

      {command: command, stdout: stdout_data, stderr_or_info: stderr_data, exit_status: exit_status }
    end


    def get_slurm_run_uuid(command, timeout)
      stderr_data_or_info = ""
      slurm_id = nil

      @ssh.open_channel do |channel|
        channel.exec(command) do |ch, success|
          unless success
            raise "Command execution failed: #{command}"
          end

          channel.on_extended_data do |_, type, data|
            stderr_data_or_info = data

            # Check for SLURM run ID pattern in stderr
            if (match = data.match(/SLURM run ID: ([a-f0-9\-]+)/))
              slurm_id = match[1]
              channel.close
            end
          end
        end
      end

      # Use a timeout to ensure the loop does not hang indefinitely
      begin
        Timeout.timeout(timeout) do
          @ssh.loop(0.1) { slurm_id.nil? }
        end
      rescue Timeout::Error
        raise "Timed out while waiting for SLURM run ID"
      end

      if slurm_id.nil?
        raise "Failed to retrieve SLURM run ID from the command output: #{stderr_data_or_info}"
      end

      slurm_id
    end

  end
end