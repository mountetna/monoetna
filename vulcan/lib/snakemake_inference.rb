class Vulcan
  class Snakemake
    module Inference
      # This class performs operations on the workspace.target_mapping, which along with the given state of the workapce,
      # helps the back-end infer which targets to build, which targets to remove, which jobs to run, etc.

      def find_buildable_targets(target_mapping, provided_params, available_files)
        # This determines the targets (files) that can be built based on the target mapping (stored in workspace.target_mapping), 
        # provided params, and available files.
        # This function really just matches the targets to the provided params and available files.
        # It is also recursive, in that it will keep adding targets to the available files, and checking again if any new targets can be built.

        # It is important to mention that this function cannot detect if files are up to date or not.
        # That means that if a file is available, it will be used as an input for a target, even if that file has changed and is outdated
        # So it is "optimistic" in this sense and just tries to build everything it can. We let snakemake handle outdated files, and it will
        # not build targets if they are outdated.
        new_target_mapping = remove_ui_targets(target_mapping)
        all_targets = Set.new
        all_files_set = available_files.to_set
        provided_params_set = provided_params.to_set
        loop do
          new_targets = match(new_target_mapping, provided_params_set, all_files_set)
          # Are there any new targets?
          break if all_targets == new_targets
          all_targets.merge(new_targets)
          all_files_set.merge(new_targets) # we add these targets to files that can be generated
        end
        all_targets.to_a
      end

      def match(target_mapping, provided_params_set, available_files_set)
        buildable = Set.new
        target_mapping.each do |target, requirements|
          # Check if all input files are available
          inputs_satisfied = requirements["inputs"].to_set.subset?(available_files_set)

          # Check if all required params are provided
          params_satisfied = requirements["params"].to_set.subset?(provided_params_set)

          # If both conditions are satisfied, add the target to buildable
          if inputs_satisfied && params_satisfied
            buildable.add(target)
          end
        end
        buildable
      end

      def find_affected_downstream_jobs(dag, scheduled_jobs)
        # It takes in a list of sorted jobs (dag), and a list of jobs that are scheduled to be run.
        # It is supposed to be able to infer which jobs downstream may be impacted by the scheduled jobs. 
        # For now, we assume a very naive approach, and we just assume affected jobs are downstream of the scheduled jobs.
        raise "No scheduled jobs provided" if scheduled_jobs.empty?
        last_matched_index = dag.rindex { |element| scheduled_jobs.include?(element) }
        raise "Cannot find any matching jobs in the dag" if last_matched_index.nil?
        # Retrieve elements after the last matched element
        dag[(last_matched_index + 1)..-1]
      end

      def ui_targets(target_mapping)
        target_mapping.select { |target, requirements| requirements["params"].include?("ui") }.keys
      end

      def file_graph(target_mapping)
        # Build reverse adjacency: input -> [outputs]
        graph = Hash.new { |h,k| h[k] = [] }
        target_mapping.each do |output, info|
          Array(info["inputs"]).each do |inp|
            graph[inp] << output
          end
        end
        graph
      end

      def downstream_nodes(graph, start_nodes)
        # Takes in an adjacency list and a list of start nodes.
        # Returns a list of all nodes that are downstream of the start nodes.
        visited = Set.new
        queue   = start_nodes.to_a.dup
      
        until queue.empty?
          node = queue.shift
          graph[node].each do |child|
            next if visited.include?(child)
            visited << child
            queue << child
          end
        end
      
        visited
      end

      def filter_ui_targets(jobs_to_run, target_mapping)
        # There is a small snakemake bug, where if a rule runs and generates a file and then
        # you quickly save a UI file - the timestamp in between the file might be too small
        # and snakemake will try to schedule a UI job. This is unlikely to ever happen in prod,
        # since by the time a config has run a UI file has been writte and run has been hit.
        # However, we just remove the UI files here to be safe.

        ui_targets = ui_targets(target_mapping)
        jobs_to_run.reject { |job| ui_targets.include?(job) }
      end

      private

      def find_targets_matching_params(target_mapping, params)
        target_mapping.select { |target, requirements| requirements["params"].any? { |param| params.include?(param) } }.keys
      end

      def remove_ui_targets(target_mapping)
        target_mapping.reject { |target, requirements| requirements["params"].include?("ui") }
      end

      module_function :find_affected_downstream_jobs, :find_buildable_targets, :match, :ui_targets, :filter_ui_targets, :find_targets_matching_params, :remove_ui_targets, :file_graph, :downstream_nodes
    end
  end
end
