require 'etna'
require 'digest'
require_relative "./vulcan_controller"
require_relative "./../../path"
require_relative './../../remote_manager'
require_relative './../../snakemake_remote_manager'
require_relative './../../snakemake_command'
require_relative './../../snakemake_inference'


class VulcanV2Controller < Vulcan::Controller

  def initialize(request, action = nil)
    super
    @remote_manager = Vulcan::RemoteManager.new(Vulcan.instance.ssh_pool)
    @snakemake_manager = Vulcan::Snakemake::RemoteManager.new(@remote_manager)
  end

  def create_workflow
    workflow = Vulcan::WorkflowV2.first(
      repo_remote_url: @params[:repo_url],
    )
    if workflow
      return success_json({'msg': "Workflow: #{workflow.name} for project: #{@params[:project_name]} already exists."})
    end
    begin
      obj = Vulcan::WorkflowV2.create(
        project_name: @params[:project_name],
        name: @params[:workflow_name],
        repo_remote_url: @params[:repo_url],
        created_at: Time.now,
        updated_at: Time.now
      )
      success_json({'workflow_id': obj.id, 'workflow_name': obj.name})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end
  # Admin
  def list_workflows
    success_json(
      workflows: Vulcan::WorkflowV2.where(
        project_name: @params[:project_name]
      ).all.map do |w|
        w.to_hash
      end
    )
  end

  def create_workspace
    workflow = Vulcan::WorkflowV2.first(id: @params[:workflow_id], project_name: @params[:project_name])
    unless workflow
      msg = "Workflow: #{@params[:workflow_name]} for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
      workspace_hash = Vulcan::Path.workspace_hash(workflow.id.to_s, @user.email)
      workspace_dir = Vulcan::Path.workspace_dir(@escaped_params[:project_name], workspace_hash)
    begin
        @remote_manager.mkdir(workspace_dir)
        @remote_manager.clone(workflow.repo_remote_url, workspace_dir)
        git_sha = @remote_manager.checkout_version(workspace_dir, @escaped_params[:git_request])
        @remote_manager.mkdir(Vulcan::Path.workspace_tmp_dir(workspace_dir))
        @remote_manager.mkdir(Vulcan::Path.workspace_output_dir(workspace_dir))
        @remote_manager.touch("#{Vulcan::Path.workspace_output_dir(workspace_dir)}/.keep")
        # @remote_manager.upload_dir(Vulcan.instance.config(:snakemake_profile_dir), workspace_dir, true) # TODO: for now we are just using the default profile
        @remote_manager.write_file(
          Vulcan::Path.dl_config(workspace_dir), 
          Vulcan::Path.dl_config_yaml(@escaped_params[:project_name], task_token, Vulcan.instance.config(:magma)[:host])
        )
        config = @remote_manager.read_yaml_file(Vulcan::Path.default_snakemake_config(workspace_dir))
        obj = Vulcan::Workspace.create(
          workflow_id: workflow.id,
          name: @params[:workspace_name],
          dag: "{#{@snakemake_manager.get_dag(workspace_dir).join(',')}}",
          target_mapping: @snakemake_manager.generate_target_mapping(workspace_dir, config),
          path: workspace_dir,
          user_email: @user.email,
          git_ref: @escaped_params[:git_request],
          git_sha: git_sha,
          created_at: Time.now,
          updated_at: Time.now
        )
        response = {
          workspace_id: obj.id,
          workflow_id: obj.workflow_id,
          vulcan_config: @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(workspace_dir)),
          dag: obj.dag
        }
        success_json(response)
      rescue => e
        @remote_manager.rmdir(workspace_dir)
        Vulcan.instance.logger.log_error(e)
        raise Etna::ServerError.new(e.message)
      end
  end

  def list_workspaces
    success_json(
      workspaces: Vulcan::Workspace.where(
        workflow_id: Vulcan::WorkflowV2.where(
          project_name: @params[:project_name]
        ).select(:id)
      ).all.map(&:to_hash)
    )
  end

  def get_workspace
    workspace = Vulcan::Workspace.first(
      id: @params[:workspace_id],
    )
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    # Fetch the last run
    last_run = Vulcan::Run.where(workspace_id: workspace.id).order(Sequel.desc(:created_at)).first
    if last_run
      job_id_hash = @snakemake_manager.parse_log_for_slurm_ids(last_run.log_path)
      slurm_status = @snakemake_manager.query_sacct(last_run.slurm_run_uuid, job_id_hash)
    end
    # Fetch the last config
    # If we dont have a config associated with the last run, just grab the most recent config
    last_config = last_run ? Vulcan::Config.where(workspace_id: workspace.id, id: last_run.config_id).first : nil
    unless last_config
      last_config = Vulcan::Config.where(workspace_id: workspace.id).order(Sequel.desc(:created_at)).first
    end

    # We update the dl_config token
    @remote_manager.write_file(
      Vulcan::Path.dl_config(workspace.path), 
      Vulcan::Path.dl_config_yaml(@escaped_params[:project_name], task_token, Vulcan.instance.config(:magma)[:host])
    )

    response = workspace.to_hash.merge({
      dag: workspace.dag,
      vulcan_config: @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(workspace.path)),
      last_config: last_config ? @remote_manager.read_json_file(last_config.path) : nil,
      last_config_id: last_config ? last_config.id : nil,
      last_run_id: last_run ? last_run.id : nil,
      last_job_status: last_run ? slurm_status : nil
    })
    success_json(response)
  end
  # Update workspace name and tags
  def update_workspace
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end

    if @params.key?(:name)
      workspace.name = @params[:name]
    end

    if @params.key?(:tags)
      formatted_tags = '{' + @params[:tags].map { |tag| "\"#{tag}\"" }.join(',') + '}'
      # Validate the format matches PostgreSQL array format
      # Our db has an old version of sequel that doesnt support the pg_array method
      safe_tags = Sequel.lit(formatted_tags).to_s
      workspace.tags = safe_tags
    end

    workspace.save
    success_json(workspace.to_hash)
  end

  def save_config
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    unless @params.key?(:uiFilesSent) && @params.key?(:paramsChanged)
      msg = "Missing required parameters: uiFilesSent and/or paramsChanged."
      raise Etna::BadRequest.new(msg)
    end
    begin
      raise Etna::TooManyRequests.new("workflow is still running...") if @snakemake_manager.snakemake_is_running?(workspace.path)
      request_params = @params[:params]
      request_params_json = to_valid_json(request_params)
      params_hash = Digest::MD5.hexdigest(request_params_json)
      config = Vulcan::Config.first(workspace_id: workspace.id, hash: params_hash)
      unless config
        # We always overwrite the default config with the incoming params
        default_config = @remote_manager.read_json_file(Vulcan::Path.default_snakemake_config(workspace.path))
        config_values = JSON.parse(default_config.to_json).merge(JSON.parse(request_params_json))
        config_path = Vulcan::Path.workspace_config_path(workspace.path, params_hash)
        @remote_manager.write_file(config_path, config_values.to_json)
        config = Vulcan::Config.create(
          workspace_id: workspace.id,
          path: config_path,
          hash: params_hash,
          created_at: Time.now,
          updated_at: Time.now
        )
      end
      @snakemake_manager.remove_existing_ui_targets(@params[:uiFilesSent], @params[:paramsChanged], workspace)
      output_files = @remote_manager.list_files(Vulcan::Path.workspace_output_dir(workspace.path)).map { |file| "output/#{file}" }
      resources_files = @remote_manager.list_files(Vulcan::Path.workspace_resources_dir(workspace.path)).map { |file| "resources/#{file}" }
      available_files = output_files + resources_files
      command = Vulcan::Snakemake::CommandBuilder.new
      command.targets = Vulcan::Snakemake::Inference.find_buildable_targets(
        workspace.target_mapping,
        request_params.keys.map(&:to_s),
        available_files
      )
      command.options = {
        config_path: config.path,
        profile_path: Vulcan::Path.profile_dir(workspace.path, "default"), # only one profile for now
        dry_run: true,
      }
      jobs_to_run = @snakemake_manager.dry_run_snakemake(workspace.path, command.build)
      jobs_to_run = Vulcan::Snakemake::Inference.filter_ui_targets(jobs_to_run, workspace.target_mapping)
      success_json(
        {
          config_id: config.id,
          scheduled: jobs_to_run,
          downstream: jobs_to_run.any? ? Vulcan::Snakemake::Inference.find_affected_downstream_jobs(workspace.dag, jobs_to_run) : [],
          params: JSON.parse(request_params_json)
        }
    )
    rescue Etna::TooManyRequests => e
      # Re-raise the TooManyRequests exception to avoid handling it as a BadRequest
      raise e
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end

  end

  def run_workflow
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    begin
      raise Etna::TooManyRequests.new("workflow is still running...") if @snakemake_manager.snakemake_is_running?(workspace.path)
      config = Vulcan::Config.where(id: @params[:config_id]).first
      unless config
        msg = "Config for workspace: #{workspace.path} does not exist."
        raise Etna::BadRequest.new(msg)
      end
      # Build snakemake command
      output_files = @remote_manager.list_files(Vulcan::Path.workspace_output_dir(workspace.path)).map { |file| "output/#{file}" }
      resources_files = @remote_manager.list_files(Vulcan::Path.workspace_resources_dir(workspace.path)).map { |file| "resources/#{file}" }
      available_files = output_files + resources_files
      params = @remote_manager.read_yaml_file(config.path).keys
      command = Vulcan::Snakemake::CommandBuilder.new
      command.targets = Vulcan::Snakemake::Inference.find_buildable_targets(workspace.target_mapping, params, available_files)
      command.options = {
        config_path: config.path,
        profile_path: Vulcan::Path.profile_dir(workspace.path, "default"), # only one profile for now
      }
      slurm_run_uuid = @snakemake_manager.run_snakemake(workspace.path, command.build)
      log = @snakemake_manager.get_snakemake_log(workspace.path, slurm_run_uuid)
      obj = Vulcan::Run.create(
        workspace_id: workspace.id,
        config_id: config.id,
        slurm_run_uuid: slurm_run_uuid,
        log_path: log,
        created_at: Time.now,
        updated_at: Time.now
      )
      success_json({run_id: obj.id})
    rescue Etna::TooManyRequests => e
      # Re-raise the TooManyRequests exception to avoid handling it as a BadRequest
      raise e
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def get_workflow_status
    workflow_run = Vulcan::Run.first(id: @params[:run_id], workspace_id: @params[:workspace_id])
    unless workflow_run
      msg = "Run for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    begin
      job_id_hash = @snakemake_manager.parse_log_for_slurm_ids(workflow_run.log_path)
      slurm_status = @snakemake_manager.query_sacct(workflow_run.slurm_run_uuid, job_id_hash)
      success_json(slurm_status)
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def cancel_workflow
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end

    # Get the specific run to cancel
    run = Vulcan::Run.first(id: @params[:run_id], workspace_id: @params[:workspace_id])
    unless run
      msg = "Run #{@params[:run_id]} not found for workspace: #{workspace.name}"
      raise Etna::BadRequest.new(msg)
    end

    begin
      # Check if the workflow is currently running
      unless @snakemake_manager.snakemake_is_running?(workspace.path)
        return success_json({
          message: "No workflow is currently running for workspace: #{workspace.name}",
          run_id: run.id
        })
      end

      # Cancel the workflow
      config = Vulcan::Config.first(id: run.config_id)
      unless config
        msg = "Config for workspace: #{workspace.path} does not exist."
        raise Etna::BadRequest.new(msg)
      end
      success = @snakemake_manager.cancel_snakemake(workspace.path, run.slurm_run_uuid)
      if success
        success_json({message: "Workflow canceled successfully", run_id: run.id})
      else
        raise Etna::BadRequest.new("Failed to cancel workflow")
      end
    rescue Etna::BadRequest => e
      raise e
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def write_files
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    raise Etna::BadRequest.new("Workspace not found") unless workspace
    files = @params[:files] || {}
    raise Etna::BadRequest.new("No files provided") if files.empty?
    begin
      files.each do |file|
        output_path = Vulcan::Path.workspace_output_dir(workspace.path)
        @remote_manager.write_file("#{output_path}#{file[:filename]}", file[:content])
      end
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
    success_json({'status': "success"})
  end

  def read_files
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    raise Etna::BadRequest.new("Workspace not found") unless workspace
    file_names = @params[:file_names] || []
    raise Etna::BadRequest.new("No files provided") if file_names.empty?
    file_contents = file_names.map do |file_name|
      output_path = Vulcan::Path.workspace_output_dir(workspace.path)
      payload = @remote_manager.read_file_to_memory("#{output_path}#{file_name}")
      utf8 = payload.dup.force_encoding('UTF-8')
      if utf8.valid_encoding?
        # safe to treat as UTF-8 text
        encoding = 'utf-8'
        content  = utf8.scrub.chomp
      else
        encoding = 'base64'
        content  = Base64.strict_encode64(payload)
      end
      {
        filename: file_name,
        content:  content,
        encoding: encoding
      }
    end

    success_json(files: file_contents)
  end

  def get_files
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    raise Etna::BadRequest.new("Workspace not found") unless workspace
    success_json({files: @remote_manager.list_files(Vulcan::Path.workspace_output_dir(workspace.path))})
  end

  def is_running
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    success_json({running: @snakemake_manager.snakemake_is_running?(workspace.path)})
  end

  def read_image
    retrieve_file(@params[:file_name], "image/png")
  end

  def download_file
    retrieve_file(@params[:file_name], "application/octet-stream", disposition: "attachment; filename=#{@params[:file_name]}")
  end

  def cluster_latency
    begin
      # Measure SSH latency using the remote_manager's measure_latency method
      median_latency = @remote_manager.measure_latency
      success_json({latency: "#{median_latency}ms"})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def delete_workspace
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace #{@params[:workspace_id]} for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end

    begin
      # Check if the workflow is currently running
      if @snakemake_manager.snakemake_is_running?(workspace.path)
        raise Etna::BadRequest.new("Cannot delete workspace while workflow is running")
      end

      # Delete all runs associated with this workspace
      Vulcan::Run.where(workspace_id: workspace.id).delete

      # Delete all configs associated with this workspace
      Vulcan::Config.where(workspace_id: workspace.id).delete

      # Delete the workspace directory on the remote server
      if @remote_manager.dir_exists?(workspace.path)
        @remote_manager.rmdir(workspace.path)
      end

      # Delete the workspace record from the database
      workspace.delete

      success_json({
        message: "Workspace #{workspace.name} deleted successfully",
        workspace_id: @params[:workspace_id]
      })
    rescue Etna::BadRequest => e
      raise e
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  private

  def retrieve_file(file_name, file_type, disposition: nil)
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    raise Etna::BadRequest.new("Workspace not found") unless workspace
    raise Etna::BadRequest.new("No file provided") if file_name.nil? || file_name.empty?
    
    output_path = Vulcan::Path.workspace_output_dir(workspace.path)
    unless @remote_manager.file_exists?("#{output_path}#{file_name}")
      raise Etna::BadRequest.new("File not found")
    end
    
    content = @remote_manager.read_file_to_memory("#{output_path}#{file_name}")
    success(content, file_type: file_type, disposition: disposition)
  end

end