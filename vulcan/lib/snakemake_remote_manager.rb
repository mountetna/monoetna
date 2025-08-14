require_relative  "./snakemake_parser"

class Vulcan
  class Snakemake
    class RemoteManager
      def initialize(remote_manager)
        @remote_manager = remote_manager
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

        sub_command = @remote_manager.build_command
          .add_raw('source ~/.bashrc') # load conda
          .add('cd', dir)
          .add('conda', 'activate', Vulcan.instance.config(:conda_env))
          .add_raw("#{snakemake_command} > #{boot_log} 2>&1 < /dev/null")

        command = @remote_manager.build_command
          .add_raw("nohup bash -lc '#{sub_command}'") 
          .background

        @remote_manager.invoke_and_close(command.to_s)
        if @remote_manager.file_exists?(boot_log, true)
          begin
            regex = /SLURM run ID: ([a-f0-9\-]+)/
            matched = nil
            3.times do |attempt| # sometimes the slurm id is not written to the boot.log immediately
              matched = @remote_manager.find_string(boot_log, regex)
              break if matched
              sleep(0.1) # 100ms delay between attempts
            end
            if matched.nil?
              file_output = @remote_manager.invoke_ssh_command("cat #{boot_log}")
              raise "Could not capture slurm id, error is: #{file_output}"
            end
            slurm_uuid = matched[1]
          rescue => e
            raise e
          end
        else
          raise "Could not create: #{boot_log}, command is: #{command.to_s}"
        end
        # We no longer need the boot.log
        @remote_manager.invoke_ssh_command(@remote_manager.build_command.add('rm', boot_log).to_s)
        slurm_uuid
      end

      def dry_run_snakemake(dir, snakemake_command)
        command = @remote_manager.build_command
          .add('conda', 'activate', Vulcan.instance.config(:conda_env))
          .add('cd', dir)
          .add_raw(snakemake_command)
        out = @remote_manager.invoke_ssh_command(command.to_s)
        get_rules_from_run(out[:stdout])
      end

      def dry_run_snakemake_files(dir, snakemake_command)
        command = @remote_manager.build_command
          .add('conda', 'activate', Vulcan.instance.config(:conda_env))
          .add('cd', dir)
          .add_raw(snakemake_command)
        out = @remote_manager.invoke_ssh_command(command.to_s)
        get_files_from_run(out[:stdout])
      end

      def get_files_from_run(snakemake_log)
        # Parse the snakemake log output to extract output file names
        # The log contains a table with columns: output_file, date, rule, log-file(s), status, plan
        # We need to extract only the output_file values where plan is "update pending"
        
        output_files = []
        
        # Split the log into lines and process each line
        snakemake_log.each_line do |line|
          line = line.strip
          next if line.empty?
          
          # Skip header lines that contain column names or separators
          next if line.match?(/^output_file\s+date\s+rule\s+log-file\(s\)\s+status\s+plan$/)
          next if line.match?(/^-+$/) # Skip separator lines
          next if line.match?(/^Building DAG of jobs\.\.\.$/) # Skip the header line
          
          # Split the line on tabs to get the proper columns
          parts = line.split("\t")
          if parts.length >= 6
            output_file = parts[0]
            plan = parts[5] # The plan column (6th column, index 5)
            
            # Only add files that have "update pending" plan
            if plan == "update pending"
              output_files << output_file unless output_file.empty?
            end
          end
        end
        
        output_files
      end

      def snakemake_is_running?(dir)
        # Snakemake drops lock files into .snakemake/locks/ when running
        # We check if the directory is empty to determine if Snakemake is running
        snakemake_locks_dir = "#{Shellwords.escape(dir)}/.snakemake/locks"
        return false unless @remote_manager.dir_exists?(snakemake_locks_dir)
        !@remote_manager.list_files(snakemake_locks_dir).empty?
      rescue => e
        Vulcan.instance.logger.error("Error checking if Snakemake is running: #{e.message}")
        false
      end

      def parse_log_for_slurm_ids(snakemake_log)
        # The slurm snakemake integration uses comments to write the name of the rule to the slurm back-end db
        # (along with the uuid). Currently --comments are disabled on the slurm instance on c4.
        # This make it very tricky to query the db and determine the mapping between rules and slurm job ids.
        # We can always query the db with the slurm snakemake uuid, and get a list of slurm jobs, but they
        # are un-identified. So, we must parse the snakemake log to figure out this mapping.
        rule_to_job_id = {}
        log = @remote_manager.read_file_to_memory(snakemake_log)
        rules = get_rules_from_run(log)
        rules.each do |rule|
          regex = /rule #{rule}:[\s\S]*?Job \d+ has been submitted with SLURM jobid (\d+)/
          match = log.match(regex)
          rule_to_job_id[rule] = match ? match[1] : nil
        end

        rule_to_job_id
      end

      def get_rules_from_run(snakemake_log)
        # TODO: refactor this

        # Extract the Job stats section up to the first double newline
        job_stats_section = snakemake_log[/Job stats:.*?(?=\n\n)/m]
        return [] unless job_stats_section
        
        # Scan for lines that have a job name followed by a count
        jobs = job_stats_section.scan(/^(\w+)\s+(\d+)$/)
        
        # Reject the 'total' entry and extract only the job names
        jobs.reject { |job, _| job.downcase == 'total' }.map { |job, _| job }
      end


      def query_sacct(slurm_uuid, job_id_hash)
        status_hash = {}
        job_id_hash.each do |job_name, job_id|
          if job_id.nil?
            status_hash[job_name] = "NOT STARTED"
          else
            command = @remote_manager.build_command
              .add('sacct', "--name=#{slurm_uuid}", '-j', job_id, '-X', '-o', 'State%20', '-n')
            out = @remote_manager.invoke_ssh_command(command.to_s)
            status_hash[job_name] = out[:stdout].empty? ? "NOT STARTED" : out[:stdout].strip
          end
        end
        status_hash
      end

      def get_snakemake_log(dir, slurm_run_uuid)
        command = @remote_manager.build_command
          .add('cd', dir)
          .add('grep', '-rl', slurm_run_uuid, Vulcan::Path::WORKSPACE_SNAKEMAKE_LOG_DIR)
        out = @remote_manager.invoke_ssh_command(command.to_s)
        log = out[:stdout].gsub("\n", "")
        "#{dir}/#{log}"
      end

      def generate_target_mapping(dir, config)
        main_snakefile = File.join(dir, "Snakefile")
        snakefiles = find_includes_in_snakefile(dir, main_snakefile)
        mapping = {}
        snakefiles.each do |file|
          file_content = @remote_manager.read_file_to_memory(file)
          parser = Vulcan::Snakemake::TargetParser.new(file_content, config)
          values = parser.parse
          mapping.merge!(values)
        end

        # Remove the all key, delete the outputs key, and remove any empty targets
        mapping.delete(:all)
        mapping
        end


      def find_includes_in_snakefile(dir, snakefile)
        included_files = [snakefile]  # Start with the main Snakefile
        result = @remote_manager.read_file_to_memory(snakefile)
        result.each_line do |line|
          # Look for 'include' statements in the Snakefile
          if line.strip.match(/^include:\s*["'](.+?)["']/)
            included_file = File.join(dir, $1)
            if @remote_manager.file_exists?(included_file)
              included_files << included_file
            else
            end
          end
        end
        included_files.uniq  # Return unique files to avoid duplicates
      end

      def get_dag(dir)
        # Snakemake REQUIRES a pipeline to have all input files EXIST, before running any of the meta commands.
        # Therefore, we adopt a convention in the config.yaml - such that all input files
        # be defined starting with the string "output/". We can then create dummy files and run
        # our meta commands.
        config_path = Vulcan::Path.default_snakemake_config(Shellwords.escape(dir))
        config = @remote_manager.read_yaml_file(config_path)

        # Identify input files and create them with "DUMMY FILE" content
        dummy_files = []
        config.each do |key, value|
          if value.is_a?(String)
            if value.start_with?('output/')
              remote_file_path = "#{Shellwords.escape(dir)}/#{Shellwords.escape(value)}"
              dummy_files.append(remote_file_path)
              @remote_manager.create_dummy_file(remote_file_path)
            end
          end
        end

        # Run snakemake --d3dag
        dag_command = @remote_manager.build_command
          .add('conda', 'activate', Vulcan.instance.config(:conda_env))
          .add('cd', dir)
          .add('snakemake', '--configfile', config_path, '--d3dag')
          .redirect_to('rulegraph.json')
        @remote_manager.invoke_ssh_command(dag_command.to_s)

        if !@remote_manager.file_exists?("#{dir}/rulegraph.json", wait: true)
          Vulcan.instance.logger.error("Failed to generate DAG: #{dag_output[:stderr_or_info]}")
          raise "Failed to generate DAG: #{dag_output[:stderr_or_info]}"
        end

        # Parse the dag and created a flattened graph
        json = @remote_manager.read_json_file("#{dir}/rulegraph.json")
        flattened_dag = d3dag_to_list(json)

        # Remove dummy files
        dummy_files.each do |value|
          @remote_manager.rm_file(value)
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

      def remove_existing_ui_targets(ui_files_written, params_changed, workspace)
        # Please see the README for more details on the halting problem.
        # TLDR; we need to delete UI step files for two cases:
        # 1. The params have changed
        # 2. New UI step files have been written

        # Only proceed if there's a reason to check (either params changed or ui files were written)
        return unless params_changed.any? || ui_files_written.any?

        # Determine what files to filter based on whether params changed or ui files were sent
        files_changed = if params_changed.any?
                                Vulcan::Snakemake::Inference.find_targets_matching_params(workspace.target_mapping, params_changed)
                              else
                                ui_files_written
                              end

        # We build a file grapho that we can only remove files that are downstream of the ui targets
        file_dag = Vulcan::Snakemake::Inference.file_graph(workspace.target_mapping)
        reachable_files = Vulcan::Snakemake::Inference.downstream_nodes(file_dag, files_changed)

        # Find the existing ui targets that are in the workspace
        existing_ui_targets = Vulcan::Snakemake::Inference.ui_targets(workspace.target_mapping)
          .select { |target| @remote_manager.file_exists?(File.join(workspace.path, target)) }

        # Remove the files
        targets_to_remove = (reachable_files & existing_ui_targets) 
        targets_to_remove.each do |target|
          Vulcan.instance.logger.info("Removing UI target: #{target}")
          @remote_manager.rm_file(File.join(workspace.path, target))
        end
      end

      def cancel_snakemake(dir, slurm_run_uuid) 
        # Find the snakemake controller process PID
        find_pid_command = @remote_manager.build_command
          .add('ps', '-eo', 'pid,cmd', '--no-headers')
          .pipe_to('awk', '/snakemake/ && /python/ && !/bash/ {print $1; exit}')
        
        pid_result = @remote_manager.invoke_ssh_command(find_pid_command.to_s)
        controller_pid = pid_result[:stdout].strip
        
        Vulcan.instance.logger.info("Found snakemake controller PID: #{controller_pid}")
        if controller_pid.empty?
          raise "No snakemake controller process found"
        end
        
        # Kill the controller process
        kill_command = @remote_manager.build_command
          .add('kill', '-2', controller_pid)
        @remote_manager.invoke_ssh_command(kill_command.to_s)
        Vulcan.instance.logger.info("Sent kill signal to PID: #{controller_pid}")
        
        # Find any left over job IDs associated with the slurm_run_uuid
        find_job_ids_command = @remote_manager.build_command
          .add('sacct', '-X', "--name=#{slurm_run_uuid}", '--format=JobID', '--noheader')
          .pipe_to('awk', 'NF')
        
        job_ids_result = @remote_manager.invoke_ssh_command(find_job_ids_command.to_s)
        job_ids = job_ids_result[:stdout].strip
        Vulcan.instance.logger.info("Found SLURM job IDs: #{job_ids}")
        
        # Cancel all the jobs using scancel
        if !job_ids.empty?
          scancel_command = @remote_manager.build_command
            .add('scancel', *job_ids.split)
          @remote_manager.invoke_ssh_command(scancel_command.to_s)
          Vulcan.instance.logger.info("Cancelled SLURM jobs for UUID: #{slurm_run_uuid}")
        else
          Vulcan.instance.logger.info("No SLURM jobs found for UUID: #{slurm_run_uuid}")
        end
      rescue => e
         Vulcan.instance.logger.error("Error canceling Snakemake workflow: #{e.message}")
      end
    end
  end
end