class Vulcan
  class WorkspaceState
    def initialize(workspace, snakemake_manager, remote_manager)
      @workspace = workspace
      @snakemake_manager = snakemake_manager
      @remote_manager = remote_manager
    end

    def state(params, config_path, available_files)
      # Find all target files that COULD be built with the given params and available files (static analysis)
      all_targets = Vulcan::Snakemake::Inference.find_buildable_targets(
        @workspace.target_mapping, 
        params, 
        available_files
      )
      # Use dry run to get the all the files, jobs that snakemake ACTUALLY plans to generate
      # This gets us the files and jobs planned and completed, also note that these are just the compute jobs and files, not the UI jobs and files.
      command = Vulcan::Snakemake::CommandBuilder.new
      command.targets = all_targets
      command.options = {
        config_path: config_path,
        profile_path: Vulcan::Path.profile_dir(@workspace.path, "default"),
        dry_run: true,
        summary: true  # Use summary to get detailed output
      }
      # === COMPUTE STATE (from snakemake dry run) ===
      # Get compute files/jobs in all three states: planned, completed, unscheduled
      # Note that these are just the compute jobs and files, not the UI jobs and files.

      dry_run_info = @snakemake_manager.dry_run_snakemake_files(@workspace.path, command.build)
      files_planned = dry_run_info[:files_scheduled].to_set.to_a
      jobs_planned = dry_run_info[:jobs_scheduled].to_set.to_a
      files_completed = dry_run_info[:files_completed].to_set.to_a
      jobs_completed = dry_run_info[:jobs_completed].to_set.to_a

      Vulcan.instance.logger.debug("Files scheduled: #{files_planned}")
      Vulcan.instance.logger.debug("Jobs scheduled: #{jobs_planned}")
      Vulcan.instance.logger.debug("Files completed: #{files_completed}")
      Vulcan.instance.logger.debug("Jobs completed: #{jobs_completed}")

      # Get all compute files and jobs from unified source (target_mapping)
      compute_files = Vulcan::Snakemake::Inference.compute_targets(@workspace.target_mapping)
      compute_job_names = Vulcan::Snakemake::Inference.compute_job_names(@workspace.target_mapping)
      
      # Calculate unscheduled compute files/jobs
      files_unscheduled = compute_files - files_completed - files_planned
      jobs_unscheduled = compute_job_names - jobs_completed - jobs_planned
      
      Vulcan.instance.logger.debug("Compute files - planned: #{files_planned}, completed: #{files_completed}, unscheduled: #{files_unscheduled}")
      Vulcan.instance.logger.debug("Compute jobs - planned: #{jobs_planned}, completed: #{jobs_completed}, unscheduled: #{jobs_unscheduled}")
      
      # === UI STATE (manually determined by file existence and staleness) ===
      # UI targets are never "planned" - only completed or unscheduled
      
      ui_state = categorize_ui_targets
      
      # Append UI state to compute state
      files_completed += ui_state[:files_completed]
      jobs_completed += ui_state[:jobs_completed]
      files_unscheduled += ui_state[:files_unscheduled]
      jobs_unscheduled += ui_state[:jobs_unscheduled]
      
      Vulcan.instance.logger.debug("Final files - planned: #{files_planned}, completed: #{files_completed}, unscheduled: #{files_unscheduled}")
      Vulcan.instance.logger.debug("Final jobs - planned: #{jobs_planned}, completed: #{jobs_completed}, unscheduled: #{jobs_unscheduled}")

      {
        available_files: available_files,
        files: {
          completed: files_completed,
          planned: files_planned,
          unscheduled: files_unscheduled
        },
        jobs: {
          completed: jobs_completed,
          planned: jobs_planned,
          unscheduled: jobs_unscheduled
        }
      }
    end

    def remove_existing_ui_targets(ui_files_written, params_changed)
      # Please see the README for more details on the halting problem.
      # TLDR; we need to delete UI step files for two cases:
      # 1. The params have changed
      # 2. New UI step files have been written

      # Only proceed if there's a reason to check (either params changed or ui files were written)
      return unless params_changed.any? || ui_files_written.any?

      # Determine what files to filter based on whether params changed or ui files were sent
      files_changed = if params_changed.any?
                              Vulcan::Snakemake::Inference.find_targets_matching_params(@workspace.target_mapping, params_changed)
                            else
                              ui_files_written
                            end

      # We build a reverse DAG that we can only remove files that are downstream of the ui targets
      file_graph = Vulcan::Snakemake::Inference.file_graph(@workspace.target_mapping)
      reachable_files = Vulcan::Snakemake::Inference.downstream_nodes(file_graph, files_changed)

      # Find the existing ui targets that are in the workspace
      existing_ui_targets = Vulcan::Snakemake::Inference.ui_targets(@workspace.target_mapping)
        .select { |target| @remote_manager.file_exists?(File.join(@workspace.path, target)) }

      # Remove the files
      targets_to_remove = (reachable_files & existing_ui_targets) 
      targets_to_remove.each do |target|
        Vulcan.instance.logger.debug("Removing UI target: #{target}")
        @remote_manager.rm_file(File.join(@workspace.path, target))
      end
    end


    def get_available_files
      output_files = @remote_manager.list_files(Vulcan::Path.workspace_output_dir(@workspace.path))
        .map { |file| "output/#{file}" }
      resources_files = @remote_manager.list_files(Vulcan::Path.workspace_resources_dir(@workspace.path))
        .map { |file| "resources/#{file}" }
      output_files + resources_files
    end

    private

    def categorize_ui_targets
      # UI targets are explicitly excluded from snakemake dry run, so we need to manually
      # check their status. A UI file is "completed" if the file exists.
      #
      # We rely on remove_existing_ui_targets to proactively delete UI files when upstream
      # params change or new files are written, so if a UI file exists, it's considered fresh.
      
      files_completed = []
      jobs_completed = []
      files_unscheduled = []
      jobs_unscheduled = []
      
      ui_targets = Vulcan::Snakemake::Inference.ui_targets(@workspace.target_mapping)
      ui_targets.each do |ui_target|
        ui_file_path = File.join(@workspace.path, ui_target)
        # Look up job name for this target (rule name extracted from target_mapping)
        job_name = Vulcan::Snakemake::Inference.rule_for_target(@workspace.target_mapping, ui_target)
        
        if @remote_manager.file_exists?(ui_file_path)
          # File exists so it is completed
          files_completed << ui_target
          jobs_completed << job_name if job_name
          Vulcan.instance.logger.debug("UI file #{ui_target} marked as completed")
        else
          # File doesn't exist so it is unscheduled
          files_unscheduled << ui_target
          jobs_unscheduled << job_name if job_name
          Vulcan.instance.logger.debug("UI file #{ui_target} marked as unscheduled: file does not exist")
        end
      end
      
      {
        files_completed: files_completed,
        jobs_completed: jobs_completed,
        files_unscheduled: files_unscheduled,
        jobs_unscheduled: jobs_unscheduled
      }
    end

  end
end 