class Vulcan
  class WorkflowV2 < Sequel::Model(:workflows)
    one_to_many :workspaces

    def to_hash
        {
          id: id,
          project: project,
          # TODO: change to name
          workflow_name: workflow_name,
          author: author,
          repo_remote_url: repo_remote_url,
          repo_local_path: repo_local_path,
          repo_tag: repo_tag,
          config: config,
          created_at: created_at,
          updated_at: updated_at
        }
    end

    def jobs
      config.keys
    end

    def valid_run_config?(run_config)
      #TODO: fix string and hash comparison
      true
    end
    #   run_config.each do |key, params|
    #     # Check if the key exists in the workflow_config
    #     unless config.key?(key)
    #       puts "Error: #{key} is not a valid top-level key in the workflow config."
    #       return false
    #     end
    #
    #     # Check if all params keys exist in the workflow_config for the given key
    #     params.each do |param_key, _|
    #       unless config[key][:params].key?(param_key)
    #         puts "Error: #{param_key} is not a valid parameter for #{key} in the workflow config."
    #         return false
    #       end
    #     end
    #   end
    #   true
    # end


  end
end
