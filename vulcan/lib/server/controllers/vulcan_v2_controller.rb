require 'etna'
require 'digest'
require_relative "./vulcan_controller"
require_relative "./../../path"
require_relative './../../remote_server_manager'
require_relative './../../snakemake'

class VulcanV2Controller < Vulcan::Controller

  def initialize(request, action = nil)
    super
    @remote_manager = Vulcan::RemoteServerManager.new(Vulcan.instance.ssh_pool)
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
        projects: "{#{@params[:projects].join(',')}}", #TODO: fix this
        name: @params[:workflow_name],
        repo_remote_url: @params[:repo_url],
        branch: @escaped_params[:branch],
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
        Sequel.lit("projects @> ARRAY['all'] OR projects @> ARRAY[?]", @params[:project_name])
      ).all.map do |w|
        w.to_hash
      end
    )
  end

  def update_workflow
  end

  def create_workspace
    workflow = Vulcan::WorkflowV2.first(id: @params[:workflow_id])
    unless workflow
      msg = "Workflow: #{@params[:workflow_name]} for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
      workspace_hash = Vulcan::Path.workspace_hash(workflow.id.to_s, @user.email)
      workspace_dir = Vulcan::Path.workspace_dir(@escaped_params[:project_name], workspace_hash)
    begin
        @remote_manager.mkdir(workspace_dir)
        @remote_manager.clone(workflow.repo_remote_url, workflow.branch, workspace_dir)
        @remote_manager.checkout_version(workspace_dir, @escaped_params[:git_version])
        @remote_manager.mkdir(Vulcan::Path.workspace_tmp_dir(workspace_dir))
        @remote_manager.upload_dir(Vulcan::Path::SNAKEMAKE_UTILS_DIR, workspace_dir, true)
        obj = Vulcan::Workspace.create(
          workflow_id: workflow.id,
          name: @params[:workspace_name],
          path: workspace_dir,
          user_email: @user.email,
          created_at: Time.now,
          updated_at: Time.now
        )
        response = {
          workspace_id: obj.id,
          workflow_id: obj.workflow_id,
          vulcan_config: @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(workspace_dir))
        }
        success_json(response)
      rescue => e
        @remote_manager.rmdir(workspace_dir)
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
  end

  def update_workspace_tag
  end

  def list_workspaces
    success_json(
      workspaces: Vulcan::Workspace.where(
        user_email: @user.email,
      ).all.map do |w|
        w.to_hash
      end
    )
  end

  def get_workspace
    workspace = Vulcan::Workspace.first(
      id: @params[:workspace_id],
      user_email: @user.email
    )
    {} unless workspace
    vulcan_config = @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(workspace.path))
    success_json(workspace.to_hash.merge(vulcan_config: vulcan_config))
  end

  def get_dag
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    begin
      success_json({dag: @remote_manager.get_dag(workspace.path)})
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
      raise Etna::TooManyRequests.new("workflow is still running...") if @remote_manager.snakemake_is_running?(workspace.path)
      config_path = Vulcan::Path.workspace_config_path(workspace.path)
      params = @params[:run] && @params[:run][:params] ? @params[:run][:params] : {} # Account for jobs with no params

      # We must always keep an updated version of the params that were used on a workflow.
      # Snakemake needs this to understand what inputs have changed. # Fetch the last params, and append new params.
      last_run = Vulcan::Run.where(workspace_id: workspace.id).order(:created_at).last
      if last_run
        old_params = @remote_manager.read_json_file(last_run.config_path)
        params = old_params.merge(params) # params take precedence over old_params
      end
      @remote_manager.write_file(config_path, params.to_json) # We must write the params of snakemake to a file to read them

      # Build snakemake commmand
      command = Vulcan::SnakemakeCommandBuilder.new
      command.options = {
        jobs: @params[:run][:jobs],
        config_path: config_path,
        profile_path: Vulcan::Path.profile_dir(workspace.path),
      }
      slurm_run_uuid = @remote_manager.run_snakemake(workspace.path, command.build)
      log = @remote_manager.get_snakemake_log(workspace.path, slurm_run_uuid)
      obj = Vulcan::Run.create(
        workspace_id: workspace.id,
        jobs: @params[:run][:jobs].join(','),
        slurm_run_uuid: slurm_run_uuid,
        config_path: config_path,
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
      job_id_hash = @remote_manager.parse_log_for_slurm_ids(workflow_run.jobs.split(','), workflow_run.log_path)
      response = @remote_manager.query_sacct(workflow_run.slurm_run_uuid, job_id_hash)
      success_json(response)
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end


  def write_files
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    raise Etna::BadRequest.new("Workspace not found") unless workspace
    raise Etna::BadRequest.new("Unsupported Content-Type: #{request.content_type}")  unless @request.content_type.start_with?('multipart/form-data')
    files = @params[:files] || []
    raise Etna:BadRequest.new("No files provided" ) if files.empty?
    begin
      files.each do |file|
        content = file[:tempfile].read()
        output_path = Vulcan::Path.workspace_output_path(workspace.path)
        @remote_manager.write_file("#{output_path}#{file[:filename]}", content)
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
    raise Etna:BadRequest.new("No files provided" ) if file_names.empty?
    file_contents = file_names.map do |file_name|
      output_path = Vulcan::Path.workspace_output_path(workspace.path)
      content = @remote_manager.read_file_to_memory("#{output_path}#{file_name}")
      {
        filename: file_name,
        content: Base64.encode64(content)
      }
    end
    success_json({files: file_contents})
  end

  def get_files
  end

  def load_params
  end

  def save_params
  end

end
