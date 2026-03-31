require 'shellwords'
require 'timeout'
require_relative "./path"
require_relative "./command_builder"

class Vulcan

  class RemoteManager
    def initialize(ssh_pool)
      @ssh_pool = ssh_pool
    end

    def build_command
      CommandBuilder.new
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
      command = build_command
        .add('mkdir', '-p', File.dirname(remote_file_path))
        .add('echo', content)
        .redirect_to(remote_file_path)

      invoke_ssh_command(command.to_s)
    end

    def read_file_to_memory(remote_file_path)
      command = build_command.add('cat', remote_file_path)
      invoke_ssh_command(command.to_s)[:stdout]
    end

    def stream_file_simple(remote_path, chunk: 256 * 1024, max_duration: 300)
      bytes_read = 0
      start_time = Time.now
      
      # Wrap the entire SFTP operation in a timeout to prevent hanging
      Timeout.timeout(max_duration) do
        @ssh_pool.with_conn do |ssh|
          sftp = ssh.sftp
          sftp.file.open(remote_path, 'r') do |file|
            loop do
              # Check timeout
              elapsed = Time.now - start_time
              if elapsed > max_duration
                Vulcan.instance.logger.warn(
                  "Stream timeout after #{elapsed.round(1)}s, #{bytes_read} bytes read from #{remote_path}"
                )
                raise Timeout::Error, "Stream exceeded #{max_duration}s limit"
              end
              
              data = file.read(chunk)
              break if data.nil? || data.empty?
              
              bytes_read += data.bytesize
              yield data
              
              # Log progress every 100MB for large files
              if bytes_read % (100 * 1024 * 1024) < chunk
                Vulcan.instance.logger.info(
                  "Streaming progress: #{(bytes_read / 1024.0 / 1024.0).round(1)}MB from #{remote_path}"
                )
              end
            end
          end
        end
      end
    rescue Timeout::Error => e
      Vulcan.instance.logger.error("Stream timeout after #{(Time.now - start_time).round(1)}s, #{bytes_read} bytes read")
      raise "Stream timeout: #{e.message}"
    rescue Net::SFTP::StatusException => e
      Vulcan.instance.logger.error("SFTP error: #{e.description}")
      raise "Remote file error: #{e.description}"
    rescue IOError, Errno::EPIPE => e
      # Client likely disconnected
      Vulcan.instance.logger.warn(
        "Client disconnected after #{bytes_read} bytes: #{e.message}"
      )
      # Don't re-raise - stream was interrupted by client
    rescue Net::SSH::Disconnect, Errno::ECONNRESET => e
      Vulcan.instance.logger.error("SSH connection error: #{e.message}")
      raise "SSH connection lost during streaming: #{e.message}"
    end

    def touch(remote_file_path)
      command = build_command.add('touch', remote_file_path)
      invoke_ssh_command(command.to_s)
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
      command = build_command.add('mkdir', '-p', dir)
      invoke_ssh_command(command.to_s)
    end

    def create_dummy_file(file_path)
      # Create a test file with the content "TEST FILE"
      command = build_command
        .add('mkdir', '-p', File.dirname(file_path))
        .add('echo', 'DUMMY FILE')
        .redirect_to(file_path)

      invoke_ssh_command(command.to_s)
    end

    def list_files(dir)
      command = build_command.add('ls', dir)
      invoke_ssh_command(command.to_s)[:stdout].split("\n")
    end

    def md5sum(file)
      command = build_command.add('md5sum', file)
      invoke_ssh_command(command.to_s)[:stdout].split("\n")
    end

    def list_dirs(dir)
      escaped_dir = dir
      # If there are no directories, return an empty array
      # "ls -d #{Shellwords.escape(dir)}/*/" throws an error if the dir is empty
      check_command = build_command
        .add('find', escaped_dir, '-maxdepth', '1', '-mindepth', '1', '-type', 'd')
        .pipe_to('wc', '-l')

      check_result = invoke_ssh_command(check_command.to_s)
      if check_result[:stdout].strip.to_i == 0
        return []
      end

      command = build_command.add('ls', '-d', "#{escaped_dir}/*/")
      result = invoke_ssh_command(command.to_s)
      result[:stdout].split("\n").map { |path| File.basename(path) }
    end

    def rmdir(dir)
      # We can solve alot of these problems by just making sure the etna-user only has permission to operate under /vulcan/app
      #TODO: handle cases: dir = /app/vulcan/workflows/ipi/
      allowed_directory = Vulcan::Path.allowed_directories.find { |allowed_dir| dir.start_with?(allowed_dir) }
      if allowed_directory.nil?
        raise ArgumentError, "Directory #{dir} cannot be deleted because it is not in the list of allowed directories."
      elsif allowed_directory == dir
        raise ArgumentError, "Cannot delete top level vulcan directory"
      else
        command = build_command.add('rm', '-r', '-f', dir)
        invoke_ssh_command(command.to_s)
      end
    end

    def rm_file(file)
      # TODO: revisit for safety
      command = build_command.add('rm', '-f', file)
      invoke_ssh_command(command.to_s)
    end

    def dir_exists?(dir)
      command = build_command
        .add('test', '-d', dir)
        .add_raw('echo "exists" || echo "not exists"')
      out = invoke_ssh_command(command.to_s)
      out[:stdout].strip == 'exists'
    end

    def clone(repo, target_dir)
      # Swap from https repo_url to form that will rely on ssh
      if repo.start_with?("https://")
        repo = repo.sub("https://github.com/", "git@github.com:") + ".git"
      end
      # For now we ignore specifically desired version, and start with the default branch
      command = build_command.add('git', 'clone', repo, target_dir)
      invoke_ssh_command(command.to_s)
    end

    def fetch_tags(dir)
      command = build_command
        .add('cd', dir)
        .add('git', 'fetch', '--tags')

      invoke_ssh_command(command.to_s)
    end

    def checkout_version(dir, sha_or_tag_or_branch)
      command = build_command
        .add('cd', dir)
        .add('git', 'checkout', sha_or_tag_or_branch)

      invoke_ssh_command(command.to_s)

      command = build_command
        .add('cd', dir)
        .add('git', 'rev-parse', '--short', 'HEAD')

      result = invoke_ssh_command(command.to_s)
      result[:stdout].chomp
    end

    def get_repo_remote_url(repo_dir)
      command = build_command
        .add('cd', repo_dir)
        .add('git', 'config', '--get', 'remote.origin.url')

      result = invoke_ssh_command(command.to_s)
      result[:stdout].chomp
    end

    def file_exists?(file_path, wait = false)
      def remote_file_exists?(file_path)
        begin
          command = build_command
            .add('test', '-f', file_path)
            .add_raw('echo "exists" || echo "not exists"')
          result = invoke_ssh_command(command.to_s)
          result[:stdout].strip == 'exists'
        rescue => e
          raise e unless e.message.include?("timed out")
        end
      end

      if wait
        max_retries = 5
        attempt = 0
        file_found = false

        while attempt <= max_retries && !file_found
          if remote_file_exists?(file_path)
            file_found = true
          else
            attempt += 1
            # Exponential backoff - sleep longer each retry
            sleep_time = 0.5 * (2 ** attempt) 
            sleep sleep_time
          end
          if file_found
            break
          end
        end
        file_found
      else
        file_found = remote_file_exists?(file_path)
      end
     
      if !file_found
        Vulcan.instance.logger.warn("File #{file_path} not found")
      end
      file_found
    end

    def file_mtime(file_path)
      # Get the modification time of a file
      return nil unless file_exists?(file_path)
      
      command = build_command
        .add('stat', '-c', '%Y', file_path) # %Y gives mtime as Unix timestamp
      result = invoke_ssh_command(command.to_s)
      
      timestamp = result[:stdout].strip.to_i
      Time.at(timestamp)
    rescue => e
      Vulcan.instance.logger.error("Error getting mtime for #{file_path}: #{e.message}")
      nil
    end

    def invoke_and_close(command)
      # Immediately run ssh command and close the ssh channel.
      @ssh_pool.with_conn do |ssh|
        ssh.open_channel do |channel|
          channel.exec(command)
          channel.close
        end
      end
    end

    def measure_latency
      # Measure SSH latency by running a simple command multiple times and taking the median
      start_time = Time.now
      command = build_command.add('echo', 'latency_test')
      result = invoke_ssh_command(command.to_s)
      end_time = Time.now
      
      latency_ms = ((end_time - start_time) * 1000).round(2)

      return latency_ms
    end
    
    def invoke_ssh_command(command, timeout: 10, retries: 5, retry_delay: 2)
      last_error  = nil
      retry_count = 0
      backoff     = retry_delay

      while retry_count < retries
        begin
          return exec_with_timeout(command, timeout)
        rescue => e
          last_error  = e
          retry_count += 1
          Vulcan.instance.logger.warn(
            "Command `#{command}` failed (#{e.class}: #{e.message}) " \
            "– attempt #{retry_count}/#{retries}, retrying in #{backoff}s…"
          )
          sleep(backoff) if retry_count < retries
          backoff *= 2
        end
      end

      raise "Command failed after #{retries} attempts. Last error: #{last_error.class}: #{last_error.message}"
    end

    private

    def exec_with_timeout(command, timeout)
      stdout_data  = ""
      stderr_data  = ""
      exit_status  = nil
      completed    = false
      timed_out    = false

      @ssh_pool.with_conn do |ssh|
        channel = ssh.open_channel do |ch|
          ch.exec(command) { |_, ok| raise "exec failed" unless ok }

          ch.on_data          { |_, data| stdout_data << data }
          ch.on_extended_data { |_, _, data| stderr_data << data }
          ch.on_request("exit-status") do |_, data|
            exit_status = data.read_long
            completed   = true
          end
          ch.on_close { completed = true }
        end

        # watcher thread only closes the one channel and marks timed_out
        watcher = Thread.new do
          sleep(timeout)
          unless completed
            channel.close
            timed_out = true
          end
        end

        until completed
          ssh.loop(0.1)
        end

        watcher.kill
        raise "Command execution timed out: #{command}" if timed_out

        if exit_status != 0
          raise "Command exited #{exit_status}: #{stderr_data}"
        end

        {
          command:        command,
          stdout:         stdout_data,
          stderr_or_info: stderr_data,
          exit_status:    exit_status
        }
      end
    end
  end
end
