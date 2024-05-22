require 'shellwords'
require 'timeout'
require_relative "./path"

class Vulcan
  class SSH
    def initialize(net_ssh_instance)
      @ssh = net_ssh_instance
    end

    def run_snakemake(dir, snakemake_command)
      # The slurm snakemake integration uses a uuid in order to track workflow submissions.
      # It logs the uuid when snakemake has properly booted up.
      # So we run the snakemake command as a background process and write its output to
      # tmp/boot.log. We then tail this file to capture this output and then immediately exit.
      # The boot.log can then me removed since snakemake is also logging this output in
      # a canonical location.
      # TODO: investigate why I must use ssh.exec! in this entire function
      # TODO: figure out how to properly manage the ssh connection

      file_to_tail = "#{dir}/#{Vulcan::Path::WORKSPACE_BOOT_LOG}"
      slurm_uuid = nil

      command = "nohup sh -c 'cd #{Shellwords.escape(dir)} && #{snakemake_command} > #{Vulcan::Path::WORKSPACE_BOOT_LOG} 2>&1' &"
      # TODO: why don't I need a  @ssh.loop here? This is the only way I can run background commands
      @ssh.exec(command)
      if file_exists(file_to_tail, true)
        regex = /SLURM run ID: ([a-f0-9\-]+)/
        begin
          slurm_uuid = tail_and_match(file_to_tail, regex, 3)
        rescue => e
          file_output = @ssh.exec!("cat #{file_to_tail}")
          require 'pry'; binding.pry
          raise "Could not capture slurm id, error is: #{file_output}"
        end
      else
        raise "Could not create: #{file_to_tail}, command is: #{command}"
      end
      # We no longer need the boot.log
      @ssh.exec!("rm #{Shellwords.escape(file_to_tail)}")
      slurm_uuid
    end

    def parse_log_for_slurm_ids(snakemake_log)
      # The slurm snakemake integration uses comments to write the name of the rule to the slurm back-end db
      # (along with the uuid). Currently --comments are disabled on the slurm instance on c4.
      # This make it very tricky to query the db and determine the mapping between rules and slurm job ids.
      # We can always query the db with the slurm snakemake uuid, and get a list of slurm jobs, but they
      # are un-identified. So, we must parse the snakemake log to figure out this mapping.
      out = read_file_to_memory(snakemake_log)
    end


    def get_snakemake_log(dir, slurm_run_uuid)
      command = "cd #{Shellwords.escape(dir)} && grep -rl #{slurm_run_uuid} #{Vulcan::Path::WORKSPACE_SNAKEMAKE_LOG_DIR}"
      out = @ssh.exec!(command)
      out.gsub("\n", "")
    end

    def write_file(remote_file_path, content)
      command = "mkdir -p $(dirname #{remote_file_path}) && echo #{Shellwords.escape(content)} > #{remote_file_path}"
      invoke_ssh_command(command)
    end

    def read_file_to_memory(remote_file_path)
      command = Shellwords.join(['cat', remote_file_path])
      result = invoke_ssh_command(command)
      result[:stdout]
    end

    def read_yaml_file(remote_file_path)
      YAML.safe_load(read_file_to_memory(remote_file_path))
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

    def file_exists(file_path, wait = False)
      def remote_file_exists?(file_path)
        result = @ssh.exec!("test -f #{file_path} && echo 'exists' || echo 'not exists'")
        result.strip == 'exists'
      end

      if wait
        max_retries = 3
        attempt = 0
        file_found = false

        while attempt <= max_retries && !file_found
          if remote_file_exists?(file_path)
            file_found = true
          else
            attempt += 1
            sleep 0.5
          end
          if file_found
            break
          end
        end
        file_found
      else
        remote_file_exists?(file_path)
      end
    end


    def tail_and_match(file_to_tail, regex, timeout_duration)
      captured_string = nil
      begin
        Timeout.timeout(timeout_duration) do
          channel_data = ""
          @ssh.exec!("tail -f #{file_to_tail}") do |channel, stream, data|
            channel_data << data
            match_data = channel_data.match(regex)
            if match_data
              captured_string = match_data[1]
              # Close the channel to stop the tail command
              channel.close
            end
          end
        end
      rescue Timeout::Error
        # Handle the timeout error if necessary
        raise "Could not tail and match. Operation timed out after #{timeout_duration} seconds."
      end
      captured_string
    end


    private

    def invoke_ssh_command(command)
      stdout_data = ""
      stderr_data = ""
      exit_status = nil
      completed = false

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
            completed = true
          end

          channel.on_close do
            completed = true
          end
        end
      end

      # Loop with a timeout to avoid long wait times
      timeout = 10 # seconds
      start_time = Time.now

      until completed || (Time.now - start_time) > timeout
        @ssh.loop(0.1) # 0.1 second loop interval
      end

      if (Time.now - start_time) > timeout
        raise "Command execution timeout: #{command}"
      end

      if exit_status != 0
        raise "Command exited with status #{exit_status}. \n Command: #{command} \n Msg: #{stderr_data}"
      end

      {command: command, stdout: stdout_data, stderr_or_info: stderr_data, exit_status: exit_status }
    end


  end
end