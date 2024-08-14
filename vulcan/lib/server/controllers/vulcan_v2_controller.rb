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

  # Admin command complete prior
  def clone_repo
    repo_name = File.basename(@escaped_params[:repo_url])
    repo_local_path = Vulcan::Path.repo_path(@escaped_params[:project_name], repo_name)
    begin
      if @remote_manager.dir_exists?("#{repo_local_path}")
        msg = "Repo: #{repo_local_path} for project: #{@params[:project_name]} already exists."
        Vulcan.instance.logger.info(msg)
        return success_json({'msg': msg})
      end
      @remote_manager.mkdir(Vulcan::Path.project_dir(@escaped_params[:project_name]))
      @remote_manager.clone(@escaped_params[:repo_url], @escaped_params[:branch], repo_local_path)
      success_json({ repo_path: repo_local_path, repo_name: repo_name})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  # Admin
  def list_repos
    begin
      project_dir = Vulcan::Path.project_dir(@escaped_params[:project_name])
      success_json(@remote_manager.list_repos_with_tags(project_dir))
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def delete_repo
    begin
      repo_path = Vulcan::Path.repo_path(@params[:project_name], @params[:repo_name])
      raise Etna::BadRequest.new("Repo #{@params[:repo_name]} does not exist") unless @remote_manager.dir_exists?(repo_path)
      @remote_manager.rmdir(repo_path)
      success_json({'msg': "Repo has been deleted"})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def pull_repo
    begin
      repo_path = Vulcan::Path.repo_path(@params[:project_name], @params[:repo_name])
      @remote_manager.fetch_tags(repo_path)
      success_json({'msg': "Repo has been pulled"})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  # User
  def publish_workflow
    workflow = Vulcan::WorkflowV2.first(
      project: @params[:project_name],
      name: @params[:workflow_name],
      repo_tag: @params[:tag]
    )
    if workflow
      return success_json({'msg': "Workflow: #{@params[:name]} for project: #{@params[:project_name]} already exists."})
    end
    tmp_dir = Vulcan::Path.tmp_dir(Vulcan::Path.tmp_hash(@escaped_params[:workflow_name], @user.email))
    begin
      # Create a temporary directory to do work inside
      @remote_manager.mkdir(tmp_dir)
      @remote_manager.clone(@escaped_params[:repo_path], @escaped_params[:branch], tmp_dir)
      @remote_manager.checkout_tag(tmp_dir, @escaped_params[:tag])
      vulcan_config = @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(tmp_dir))
      # TODO: maybe run a validation on the snakefile, config and the vulcan_config
      obj = Vulcan::WorkflowV2.create(
        project: @escaped_params[:project_name],
        name: @escaped_params[:workflow_name],
        author: @escaped_params[:author],
        repo_remote_url: @remote_manager.get_repo_remote_url(@escaped_params[:repo_path]),
        repo_path: @escaped_params[:repo_path],
        branch: @escaped_params[:branch],
        vulcan_config: JSON.pretty_generate(vulcan_config),
        repo_tag: @escaped_params[:tag],
        created_at: Time.now,
        updated_at: Time.now
      )
      @remote_manager.rmdir(tmp_dir)
      success_json({'workflow_id': obj.id, 'workflow_name': obj.name})
    rescue => e
      @remote_manager.rmdir(tmp_dir)
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def list_workflows
    success_json(
      workflows: Vulcan::WorkflowV2.where(
        project: @params[:project_name]
      ).all.map do |w|
        w.to_hash
      end
    )
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
        @remote_manager.clone(workflow.repo_path, workflow.branch, workspace_dir)
        @remote_manager.checkout_tag(workspace_dir, workflow.repo_tag)
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
          vulcan_config: workflow.vulcan_config
        }
        success_json(response)
      rescue => e
        @remote_manager.rmdir(workspace_dir)
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
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
    success_json(
      workspace: Vulcan::Workspace.first(
        id: @params[:workspace_id],
        user_email: @user.email
      ).to_hash_with_vulcan_config)
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
      config_path = Vulcan::Path.workspace_config_path(workspace.path)
      @remote_manager.write_file(config_path, @params[:run][:params].to_json) # we must write the params of snakemake to a file to read them
      command = Vulcan::SnakemakeCommandBuilder.new
      command.options = {
        run_until: @params[:run][:jobs].last,
        config_path: config_path,
        profile_path: Vulcan::Path.profile_dir(workspace.path),
      }
      slurm_run_uuid = @remote_manager.run_snakemake(workspace.path, command.build)
      log = @remote_manager.get_snakemake_log(workspace.path, slurm_run_uuid)
      obj = Vulcan::Run.create(
        workspace_id: workspace.id,
        slurm_run_uuid: slurm_run_uuid,
        config_path: config_path,
        log_path: log,
        created_at: Time.now,
        updated_at: Time.now
      )
      success_json({run_id: obj.id})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def get_workflow_status
    workflow_run = Vulcan::Run.first(id: @params[:run_id], workspace_id: @params[:workspace_id])
    if workflow_run
      response = {}
      begin
        config = @remote_manager.read_json_file(workflow_run.config_path)
        job_id_hash = @remote_manager.parse_log_for_slurm_ids(config.keys, workflow_run.snakemake_log_path)
        response = @remote_manager.query_sacct(workflow_run.slurm_run_uuid, job_id_hash)
      rescue => e
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
    else
      msg = "Workflow for project: #{@params[:project_name]} does not exist."
      response = {'Warning': msg}
    end
    success_json(response)
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

end
