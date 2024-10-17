class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow_v2, :key => :workflow_id
    one_to_many :configs

    def to_hash
      {
        workspace_id: id,
        workspace_name: name,
        workflow_id: workflow_id,
        workflow_name: workflow.name,
        user_email: user_email,
        tags: tags,
        workspace_path: path,
        created_at: created_at,
        updated_at: updated_at
      }
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


    def find_buildable_targets(provided_params, available_files)
      all_targets = Set.new
      all_files_set = available_files.to_set
      provided_params_set = provided_params.to_set
      loop do
        new_targets = match(target_mapping, provided_params_set, all_files_set)
        # Are there any new targets?
        break if all_targets == new_targets
        all_targets.merge(new_targets)
        all_files_set.merge(new_targets) # we add these targets to files that can be generated
      end
      all_targets.to_a
    end
  end
end
