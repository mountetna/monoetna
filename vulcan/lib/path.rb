# frozen_string_literal: true

class Vulcan
  
  module Path
    WORKSPACE_BOOT_LOG = "tmp/boot.log"
    WORKSPACE_SNAKEMAKE_LOG_DIR = ".snakemake/log/"

    # Root paths

    def self.base_dir
      Vulcan.instance.config(:base_dir)
    end
    
    def self.workspace_base_dir
      "#{base_dir}/workspace"
    end

    # Workspace paths

    def self.workspace_tmp_dir(workspace_dir)
      "#{workspace_dir}/tmp/" 
    end

    def self.workspace_dir(project_name, hash)
      "#{self.workspace_base_dir}/#{project_name}/#{hash}"
    end

    def self.workspace_hash(workflow_id, user_email)
      Digest::MD5.hexdigest(workflow_id + Time.now.to_s + user_email)
    end

    def self.vulcan_config(workspace_dir)
      "#{workspace_dir}/vulcan_config.yaml"
    end

    def self.dl_config(workspace_dir)
      "#{workspace_dir}/dl_config.yaml"
    end

    def self.default_snakemake_config(workspace_dir)
      "#{workspace_dir}/default-config.json"
    end

    def self.workspace_config_dir(workspace_dir)
      "#{workspace_dir}/configs/"
    end

    def self.workspace_config_path(workspace_dir, md5_hash)
      self.workspace_config_dir(workspace_dir) + "#{md5_hash}.json"
    end

    def self.workspace_output_dir(workspace_dir)
      "#{workspace_dir}/output/"
    end

    def self.workspace_resources_dir(workspace_dir)
      "#{workspace_dir}/resources/"
    end

    def self.profile_dir(workspace_dir, profile_name)
      "#{workspace_dir}/profiles/#{profile_name}/"
    end

    # Other
    def self.dl_config_yaml(project_name, token, magma_url)
      yaml_content = <<~YAML
      project_name: "#{project_name}"
      token: "#{token}"
      magma_url: "#{magma_url}"
    YAML
    end

    def self.allowed_directories
      [self.workspace_base_dir]
    end

  end
end
