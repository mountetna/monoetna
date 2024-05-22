# frozen_string_literal: true

class Vulcan
  module Path
    # Potentially move these to .env
    WORKFLOW_BASE_DIR = "/app/vulcan/workflows"
    WORKSPACE_BASE_DIR = "/app/vulcan/workspace"
    WORKSPACE_TMP_DIR = "tmp/"
    WORKSPACE_BOOT_LOG = "tmp/boot.log"
    WORKSPACE_SNAKEMAKE_LOG_DIR = ".snakemake/log/"
    WORKSPACE_PROFILE_DIR = "snakemake_utils/profiles/test/"
    TMPDIR = "/tmp/vulcan"
    ALLOWED_DIRECTORIES = [WORKSPACE_BASE_DIR, WORKFLOW_BASE_DIR, TMPDIR]

    # TODO: revisit
    SNAKEMAKE_UTILS_DIR = "/app/snakemake_utils" # This is local

    def self.project_dir(project_name)
      "#{WORKFLOW_BASE_DIR}/#{project_name}"
    end

    def self.repo_local_path(project_name, repository_name)
      "#{project_dir(project_name)}/#{repository_name}"
    end

    def self.workspace_dir(project_name, hash)
      "#{WORKSPACE_BASE_DIR}/#{project_name}/#{hash}"
    end

    def self.workspace_hash(workflow_id, user_email)
      Digest::MD5.hexdigest(workflow_id + Time.now.to_s + user_email)
    end

    def self.workspace_tmp_dir(workspace)
      "#{workspace}/#{WORKSPACE_TMP_DIR}" # create a tmp directory
    end

    def self.tmp_hash(workflow_name, user_email)
      Digest::MD5.hexdigest(workflow_name + Time.now.to_s + user_email)
    end

    def self.tmp_dir(tmp_hash)
      "#{TMPDIR}/#{tmp_hash}"
    end

    def self.metis_mirror_path(workspace_dir)
      "#{workspace_dir}/metis_output/"
    end

    def self.workspace_config_path(workspace_path)
      "#{workspace_path}/run_config/#{Time.now.to_i}.json"
    end

    def self.profile_dir(workspace_path)
      "#{workspace_path}/#{WORKSPACE_PROFILE_DIR}"
    end

  end
end
