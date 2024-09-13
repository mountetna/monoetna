require 'shellwords'
require 'timeout'
require_relative "./path"

class Vulcan

  class RemoteServerManager
    def initialize(ssh_pool)
      @ssh_pool = ssh_pool
    end

    def run_snakemake(dir, snakemake_command)
      # The slurm snakemake integration uses a uuid in order to track workflow submissions.
      # It logs the uuid when snakemake has properly booted up.
      # So we run the snakemake command as a background process and write its output to
      # tmp/boot.log. We then tail this file to capture this output and then immediately exit.
      # The boot.log can then me removed since snakemake is also logging this output in
      # a canonical location.

      boot_log = "#{dir}/#{Vulcan::Path::WORKSPACE_BOOT_LOG}"
      slurm_uuid = nil

      command = "cd #{Shellwords.escape(dir)} && #{snakemake_command} > #{Vulcan::Path::WORKSPACE_BOOT_LOG}"
      invoke_background_ssh_command(command)
      if file_exists?(boot_log, true)
        begin
          regex = /SLURM run ID: ([a-f0-9\-]+)/
          matched = find_string(boot_log, regex)
          if matched.nil?
            file_output = invoke_ssh_command("cat #{boot_log}")
            raise "Could not capture slurm id, error is: #{file_output}"
          end
          slurm_uuid = matched[1]
        rescue => e
          raise e
        end
      else
        raise "Could not create: #{boot_log}, command is: #{command}"
      end
      # We no longer need the boot.log
      invoke_ssh_command("rm #{Shellwords.escape(boot_log)}")
      slurm_uuid
    end


    def snakemake_is_running?(dir)
      # Snakemake uses a lock file to ensure that only one instance of Snakemake is running in a particular working directory at a time.
      # You will run into lock errors if you try to invoke snakemake before the last process has finished running.
      # The function checks if there are any open files in the .snakemake/log directory.
      # If any files are open, it assumes that a Snakemake process is running.
      snakemake_log_dir = "#{Shellwords.escape(dir)}/.snakemake/log"
      return false unless dir_exists?(snakemake_log_dir)
      command = "lsof +D #{snakemake_log_dir}"
      begin
        result = invoke_ssh_command(command)
        # Return true if files are open, false otherwise
        result[:exit_status] == 0
      rescue RuntimeError => e
        # invoke_ssh_command raises an error for any return status != 0,
        # and for some reason lsof seems to always return 1 even when there
        # are open files. When there are files open lsof displays:
        # COMMAND PID USER   FD   TYPE DEVICE SIZE/OFF NODE NAME
        # so we just search for this string in stdout.
        return true if e.message.include?("COMMAND")
        return false if e.message.include?("Command exited with status 1")
        raise "Error checking if Snakemake is running: #{e.message}"
      end
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

    def parse_log_for_slurm_ids(rules, snakemake_log)
      # The slurm snakemake integration uses comments to write the name of the rule to the slurm back-end db
      # (along with the uuid). Currently --comments are disabled on the slurm instance on c4.
      # This make it very tricky to query the db and determine the mapping between rules and slurm job ids.
      # We can always query the db with the slurm snakemake uuid, and get a list of slurm jobs, but they
      # are un-identified. So, we must parse the snakemake log to figure out this mapping.
      rule_to_job_id = {}
      log = read_file_to_memory(snakemake_log)

      rules.each do |rule|
        regex = /rule #{rule}:[\s\S]*?Job \d+ has been submitted with SLURM jobid (\d+)/
        match = log.match(regex)
        rule_to_job_id[rule] = match ? match[1] : nil
      end

      rule_to_job_id
    end

    def query_sacct(slurm_uuid, job_id_hash)
      status_hash = {}
      job_id_hash.each do |job_name, job_id|
        if job_id.nil?
          status_hash[job_name] = "NOT STARTED"
        else
          command = "sacct --name=#{Shellwords.escape(slurm_uuid)} -j #{Shellwords.escape(job_id)} -X -o State -n"
          out = invoke_ssh_command(command)
          status_hash[job_name] = out[:stdout].empty? ?  "NOT STARTED" : out[:stdout].strip
        end
      end
      status_hash
    end


    def get_snakemake_log(dir, slurm_run_uuid)
      command = "cd #{Shellwords.escape(dir)} && grep -rl #{slurm_run_uuid} #{Vulcan::Path::WORKSPACE_SNAKEMAKE_LOG_DIR}"
      out = invoke_ssh_command(command)
      log = out[:stdout].gsub("\n", "")
      "#{dir}/#{log}"
    end

    def get_dag(dir)
      # Snakemake REQUIRES a pipeline to have all input files EXIST, before running any of the meta commands.
      # Therefore, we adopt a convention in the config.yaml - such that all input files
      # be defined starting with the string "output/". We can then create dummy files and run
      # our meta commands.

      config_path = "#{Shellwords.escape(dir)}/config.yaml"
      config = read_yaml_file(config_path)

      # Identify input files and create them with "DUMMY FILE" content
      dummy_files = []
      config.each do |key, value|
        if value.is_a?(String)
          if value.start_with?('output/')
          remote_file_path = "#{Shellwords.escape(dir)}/#{Shellwords.escape(value)}"
          dummy_files.append(remote_file_path)
          create_dummy_file(remote_file_path)
          end
        end
      end

      # Run snakemake --d3dag
      dag_command = "cd #{Shellwords.escape(dir)} && snakemake --d3dag > rulegraph.json"
      dag_output = invoke_ssh_command(dag_command)
      if dag_output[:exit_status] != 0
        raise "Failed to generate DAG: #{dag_output[:stderr_or_info]}"
      end

      # Parse the dag and created a flattened graph
      json = read_json_file("#{dir}/rulegraph.json")
      flattened_dag = d3dag_to_list(json)

      # Remove dummy files
      dummy_files.each do |value|
        rm_file(value)
      end
      flattened_dag
    end

    def d3dag_to_list(json)
      # Map job IDs to rule labels
      id_to_label = json['nodes'].map { |node| [node['id'], node['value']['label']] }.to_h

      # Create dependencies hash
      dependencies = Hash.new { |hash, key| hash[key] = [] }
      json['links'].each do |link|
        dependencies[id_to_label[link['v']]] << id_to_label[link['u']]
      end

      # Topological sort
      sorted_rules, visited = [], {}
      visit = lambda do |rule|
        next if visited[rule]
        visited[rule] = true
        dependencies[rule].each { |dep| visit.call(dep) }
        sorted_rules << rule
      end

      id_to_label.values.each { |rule| visit.call(rule) unless visited[rule] }
      sorted_rules.reject { |rule| rule == "all" }
    end

    def create_dummy_file(file_path)
      # Create a test file with the content "TEST FILE"
      command = "mkdir -p $(dirname #{file_path}) && echo 'DUMMY FILE' > #{file_path}"
      invoke_ssh_command(command)
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

    def list_repos_with_tags(dir)
      repos = list_dirs(dir)
      repo_tag_map = {}
      repos.each do |repo|
        path = "#{dir}/#{repo}"
        command = "cd #{Shellwords.escape(path)} && git tag"
        out = invoke_ssh_command(command)
        repo_tag_map[repo] = out[:stdout].split("\n")
      end
      repo_tag_map
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
      # TODO: revisit
      command = Shellwords.join(["rm","-f", file])
      invoke_ssh_command(command)
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

    private

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