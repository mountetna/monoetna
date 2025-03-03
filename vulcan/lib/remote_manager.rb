require 'shellwords'
require 'timeout'
require_relative "./path"

class Vulcan

  class RemoteManager
    def initialize(ssh_pool)
      @ssh_pool = ssh_pool
    end

    def find_string(file, regex, attempts = 3)
      found = nil
      attempts.times do
        matched = read_file_to_memory(file).match(regex)
        if matched
          found = matched
        else
          sleep 0.5
        end
      end
      found
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
      YAML.safe_load(read_file_to_memory(remote_file_path), permitted_classes: [Symbol])
    end

    def read_json_file(remote_file_path)
      JSON.parse(read_file_to_memory(remote_file_path))
  end

    def upload_dir(local_dir, remote_dir, recurse)
      @ssh_pool.with_conn do |ssh|
        ssh.scp.upload!(local_dir, remote_dir, recursive: recurse)
      end
    end

    def mkdir(dir)
      # Make project directory if it doesnt exist
      command = Shellwords.join(["mkdir", "-p", dir])
      invoke_ssh_command(command)
    end

    def create_dummy_file(file_path)
      # Create a test file with the content "TEST FILE"
      command = "mkdir -p $(dirname #{file_path}) && echo 'DUMMY FILE' > #{file_path}"
      invoke_ssh_command(command)
    end

    def list_files(dir, append_dir: false)
      command = "ls #{Shellwords.escape(dir)}"
      out = invoke_ssh_command(command)
      out[:stdout].split("\n")
    end

    def md5sum(file)
      command = "md5sum #{Shellwords.escape(file)}"
      out = invoke_ssh_command(command)
      out[:stdout].split("\n")
    end

    def list_dirs(dir)
      escaped_dir = Shellwords.escape(dir)
      # If there are no directories, return an empty array
      # "ls -d #{Shellwords.escape(dir)}/*/" throws an error if the dir is empty
      check_dirs_command = "find #{escaped_dir} -maxdepth 1 -mindepth 1 -type d | wc -l"
      check_result = invoke_ssh_command(check_dirs_command)
      if check_result[:stdout].strip.to_i == 0
        return []
      end
      command = "ls -d #{Shellwords.escape(dir)}/*/"
      result = invoke_ssh_command(command)
      result[:stdout].split("\n").map { |path| File.basename(path) }
    end

    def rmdir(dir)
      # We can solve alot of these problems by just making sure the etna-user only has permission to operate under /vulcan/app
      #TODO: handle cases: dir = /app/vulcan/workflows/ipi/
      allowed_directory = Vulcan::Path::ALLOWED_DIRECTORIES.find { |allowed_dir| dir.start_with?(allowed_dir) }
      if allowed_directory.nil?
        raise ArgumentError, "Directory #{dir} cannot be deleted because it is not in the list of allowed directories."
      elsif allowed_directory == dir
        raise ArgumentError, "Cannot delete top level vulcan directory"
      else
        command = Shellwords.join(["rm", "-r", "-f", dir])
        invoke_ssh_command(command)
      end
    end

    def rm_file(file)
      # TODO: revisit for safety
        command = Shellwords.join(["rm","-f", file])
      invoke_ssh_command(command)
    end

    def dir_exists?(dir)
      command = "[ -d #{dir} ] && echo 'Directory exists.' || echo 'Directory does not exist.'"
      out = invoke_ssh_command(command)
      out[:stdout].chomp == 'Directory exists.'
    end

    def clone(repo, target_dir)
      # For now we ignore branch and assume the default branch
      command = Shellwords.join(['git', 'clone', repo, target_dir])
      invoke_ssh_command(command)
    end

    def fetch_tags(dir)
      command = "cd #{Shellwords.escape(dir)} && git fetch --tags"
      invoke_ssh_command(command)
    end

    def checkout_version(dir, sha_or_tag)
      command = "cd #{Shellwords.escape(dir)} && git checkout #{Shellwords.escape(sha_or_tag)}"
      invoke_ssh_command(command)
    end

    def get_repo_remote_url(repo_dir)
      command = "cd #{Shellwords.escape(repo_dir)} && git config --get remote.origin.url"
      result = invoke_ssh_command(command)
      result[:stdout].chomp
    end

    def file_exists?(file_path, wait = false)
      def remote_file_exists?(file_path)
        begin
          command = "test -f #{file_path} && echo 'exists' || echo 'not exists'"
          result = invoke_ssh_command(command, 1)
          result[:stdout].strip == 'exists'
        rescue => e
          raise e unless e.message.include?("timed out")
        end
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

    def invoke_background_ssh_command(command)
      # Runs a ssh command as a async background process.
      # Combines standard err and standard out to the same stream
      # Immediately closes the ssh channel.
      wrapped_command = "nohup sh -c '#{command} 2>&1' &"
      @ssh_pool.with_conn do |ssh|
        ssh.open_channel do |channel|
          channel.exec(wrapped_command)
          channel.close
          end
        end
    end

    def invoke_ssh_command(command, timeout = 10)
      # This function runs a async command and keeps polling until the command has completed
      # or until the timeout occurs.
      # This function gathers metadata about a command, so is useful when you want
      # explicit return information.
      stdout_data = ""
      stderr_data = ""
      exit_status = nil
      completed = false

      @ssh_pool.with_conn do |ssh|
        ssh.open_channel do |channel|
          channel.exec(command) do |ch, success|
            unless success
              raise "Command execution failed: #{command}"
            end

            ch.on_data do |_, data|
              stdout_data += data
            end

            ch.on_extended_data do |_, _, data|
              stderr_data += data
            end

            ch.on_request("exit-status") do |_, data|
              exit_status = data.read_long
              completed = true
            end

            ch.on_close do
              completed = true
            end
          end
        end

        # Start a separate thread to monitor the timeout
        timeout_thread = Thread.new do
          sleep timeout
          unless completed
            ssh.close
            raise "Command execution timed out: #{command}"
          end
        end

        until completed
          ssh.loop(0.1) # 0.1 second loop interval
        end

        # Ensure the timeout thread is terminated
        timeout_thread.kill

        if exit_status != 0
          raise "Command exited with status #{exit_status}. \n Command: #{command} \n Msg: #{stderr_data} \n Stdout: #{stdout_data}"
        end
      end
      {command: command, stdout: stdout_data, stderr_or_info: stderr_data, exit_status: exit_status }
    end

  end
end