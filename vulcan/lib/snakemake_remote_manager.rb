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

        command = "cd #{Shellwords.escape(dir)} && #{snakemake_command} > #{Vulcan::Path::WORKSPACE_BOOT_LOG}"
        @remote_manager.invoke_background_ssh_command(command)
        if @remote_manager.file_exists?(boot_log, true)
          begin
            regex = /SLURM run ID: ([a-f0-9\-]+)/
            matched = @remote_manager.find_string(boot_log, regex)
            if matched.nil?
              file_output = @remote_manager.invoke_ssh_command("cat #{boot_log}")
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
        @remote_manager.invoke_ssh_command("rm #{Shellwords.escape(boot_log)}")
        slurm_uuid
      end

      def dry_run_snakemake(dir, snakemake_command)
        command = "cd #{Shellwords.escape(dir)} && #{snakemake_command}"
        out = @remote_manager.invoke_ssh_command(command)
        get_rules_from_run(out[:stdout])
      end

      def snakemake_is_running?(dir)
        # Snakemake uses a lock file to ensure that only one instance of Snakemake is running in a particular working directory at a time.
        # You will run into lock errors if you try to @remote_manager.invoke snakemake before the last process has finished running.
        # The function checks if there are any open files in the .snakemake/log directory.
        # If any files are open, it assumes that a Snakemake process is running.
        snakemake_log_dir = "#{Shellwords.escape(dir)}/.snakemake/log"
        return false unless @remote_manager.dir_exists?(snakemake_log_dir)
        command = "lsof +D #{snakemake_log_dir}"
        begin
          result = @remote_manager.invoke_ssh_command(command)
          # Return true if files are open, false otherwise
          result[:exit_status] == 0
        rescue RuntimeError => e
          # @remote_manager.invoke_ssh_command raises an error for any return status != 0,
          # and for some reason lsof seems to always return 1 even when there
          # are open files. When there are files open lsof displays:
          # COMMAND PID USER   FD   TYPE DEVICE SIZE/OFF NODE NAME
          # so we just search for this string in stdout.
          return true if e.message.include?("COMMAND")
          return false if e.message.include?("Command exited with status 1")
          raise "Error checking if Snakemake is running: #{e.message}"
        end
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
            command = "sacct --name=#{Shellwords.escape(slurm_uuid)} -j #{Shellwords.escape(job_id)} -X -o State -n"
            out = @remote_manager.invoke_ssh_command(command)
            status_hash[job_name] = out[:stdout].empty? ? "NOT STARTED" : out[:stdout].strip
          end
        end
        status_hash
      end

      def get_snakemake_log(dir, slurm_run_uuid)
        command = "cd #{Shellwords.escape(dir)} && grep -rl #{slurm_run_uuid} #{Vulcan::Path::WORKSPACE_SNAKEMAKE_LOG_DIR}"
        out = @remote_manager.invoke_ssh_command(command)
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
        dag_command = "cd #{Shellwords.escape(dir)} && snakemake --configfile #{config_path} --d3dag > rulegraph.json"
        dag_output = @remote_manager.invoke_ssh_command(dag_command)
        if dag_output[:exit_status] != 0
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
    end
  end
end