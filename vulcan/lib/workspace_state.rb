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
      command = Vulcan::Snakemake::CommandBuilder.new
      command.targets = all_targets
      command.options = {
        config_path: config_path,
        profile_path: Vulcan::Path.profile_dir(@workspace.path, "default"),
        dry_run: true,
        summary: true  # Use summary to get detailed output
      }
      dry_run_info = @snakemake_manager.dry_run_snakemake_files(@workspace.path, command.build)
      files_planned = dry_run_info[:files_scheduled].to_set.to_a
      jobs_planned = dry_run_info[:jobs_scheduled].to_set.to_a
      files_completed = dry_run_info[:files_completed].to_set.to_a
      jobs_completed = dry_run_info[:jobs_completed].to_set.to_a

      Vulcan.instance.logger.debug("Files scheduled: #{files_planned}")
      Vulcan.instance.logger.debug("Jobs scheduled: #{jobs_planned}")
      Vulcan.instance.logger.debug("Files completed: #{files_completed}")
      Vulcan.instance.logger.debug("Jobs completed: #{jobs_completed}")

      file_graph = Vulcan::Snakemake::Inference.file_graph(@workspace.target_mapping)
      files_unscheduled = file_graph.keys - files_completed - files_planned
      jobs_unscheduled = @workspace.dag.keys - jobs_completed - jobs_planned

      Vulcan.instance.logger.debug("Unscheduled downstream files: #{files_unscheduled}")
      Vulcan.instance.logger.debug("Unscheduled downstream jobs: #{jobs_unscheduled}")

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

  end
end 