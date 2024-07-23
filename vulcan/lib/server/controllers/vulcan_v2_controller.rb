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
    repo_local_path = Vulcan::Path.repo_local_path(@escaped_params[:project_name], repo_name)
    begin
      if @remote_manager.dir_exists?("#{repo_local_path}")
        msg = "Repo: #{repo_local_path} for project: #{@params[:project_name]} already exists."
        Vulcan.instance.logger.info(msg)
        return success_json({'Warning': msg})
      end
      @remote_manager.mkdir(Vulcan::Path.project_dir(@escaped_params[:project_name]))
      @remote_manager.clone(@escaped_params[:repo_url], @escaped_params[:branch], repo_local_path)
      success_json({repo_local_path: repo_local_path, repo_name: repo_name})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  # Admin
  def list_repos
    begin
      project_dir = Vulcan::Path.project_dir(@escaped_params[:project_name])
      success_json({dirs: @remote_manager.list_dirs(project_dir)})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def delete_repo
    require 'pry' ; binding.pry
  end

  # User
  def publish_workflow
    workflow = Vulcan::WorkflowV2.first(
      project: @params[:project_name],
      workflow_name: @params[:workflow_name],
      repo_tag: @params[:tag]
    )
    if workflow
      msg = "Workflow: #{@params[:workflow_name]} for project: #{@params[:project_name]} already exists."
      response = {'Warning': msg}
    else
      begin
        repo_name = File.basename(@escaped_params[:repo_local_path])
        # Create a temporary directory to do work inside
        tmp_dir = Vulcan::Path.tmp_dir(Vulcan::Path.tmp_hash(@escaped_params[:workflow_name], @user.email))
        @remote_manager.mkdir(tmp_dir)
        @remote_manager.clone(@escaped_params[:repo_local_path], @escaped_params[:branch], tmp_dir)
        @remote_manager.checkout_tag(tmp_dir, @escaped_params[:tag])
        config = @remote_manager.read_yaml_file("#{tmp_dir}/vulcan_config.yaml")
        # TODO: run a validation on the snakefile, config and the vulcan_config
        obj = Vulcan::WorkflowV2.create(
          project: @escaped_params[:project_name],
          workflow_name: @escaped_params[:workflow_name],
          author: @escaped_params[:author], #TODO: remove author
          repo_remote_url: @remote_manager.get_repo_remote_url(@escaped_params[:repo_local_path]),
          repo_local_path: @escaped_params[:repo_local_path],
          repo_tag: @escaped_params[:tag],
          config: config.to_json,
          created_at: Time.now,
          updated_at: Time.now
        )
        @remote_manager.rmdir(tmp_dir, Vulcan::Path::ALLOWED_DIRECTORIES)
        response = {'workflow_id': obj.id, 'workflow_name': obj.workflow_name}
      rescue => e
        @remote_manager.rmdir(tmp_dir, Vulcan::Path::ALLOWED_DIRECTORIES)
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
    end
    success_json(response)
  end

  def update_workflow
    success_json({'it works!': true})
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
    response = {}
    if workflow
      workspace_hash = Vulcan::Path.workspace_hash(workflow.id.to_s, @user.email)
      workspace_dir = Vulcan::Path.workspace_dir(@escaped_params[:project_name], workspace_hash)
      begin
        @remote_manager.mkdir(workspace_dir)
        # TODO: we probably want to store the name of the "master" branch
        @remote_manager.clone(workflow.repo_local_path, "main", workspace_dir)
        @remote_manager.checkout_tag(workspace_dir, workflow.repo_tag)
        @remote_manager.mkdir("#{workspace_dir}/tmp/") # create a tmp directory
        #@remote_manager.mkdir(metis_mirror_path(workspace_dir)) # location of output files that will get mirrored back to metis
        @remote_manager.upload_dir(Vulcan::Path::SNAKEMAKE_UTILS_DIR, workspace_dir, true)
        obj = Vulcan::Workspace.create(
          workflow_id: workflow.id,
          path: workspace_dir,
          user_email: @user.email,
          created_at: Time.now,
          updated_at: Time.now
        )
        response = {
          workspace_id: obj.id,
          workflow_id: obj.workflow_id,
          workflow_config: workflow.config
        }
      rescue => e
        @remote_manager.rmdir(workspace_dir, Vulcan::Path::ALLOWED_DIRECTORIES)
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
    else
      msg = "Workflow: #{@params[:workflow_name]} for project: #{@params[:project_name]} does not exist."
      response = {'Warning': msg}
    end
    success_json(response)
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
    # Include last run
    success_json(
      workspace: Vulcan::Workspace.first(
        id: @params[:workspace_id],
        user_email: @user.email
      ).to_hash)
  end

  def run_workflow
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    if workspace
      response = {}
      begin
        run_config = @params[:run]
        if workspace.workflow_v2.valid_run_config?(run_config)
          config_path = Vulcan::Path.workspace_config_path(workspace.path)
          @remote_manager.write_file(config_path, run_config.to_json) # we must write the params of snakemake to a file to read them
          command = Vulcan::SnakemakeCommandBuilder.new
          command.options = {
            run_until: run_config.keys.last,
            config_path: config_path,
            profile_path: Vulcan::Path.profile_dir(workspace.path),
          }
          slurm_run_uuid = @remote_manager.run_snakemake(workspace.path, command.build)
          log = @remote_manager.get_snakemake_log(workspace.path, slurm_run_uuid)
          obj = Vulcan::Run.create(
            workspace_id: workspace.id,
            slurm_run_uuid: slurm_run_uuid,
            config_path: config_path[config_path.index('run_config')..-1],
            log_path: log,
            created_at: Time.now,
            updated_at: Time.now
          )
          response = {run_id: obj.id}
        else
          msg = "Workflow run config is not valid"
          response = {'Error': msg}
        end
      rescue => e
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
    else
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      response = {'Warning': msg}
    end
    success_json(response)
  end

  def get_workflow_status
    workflow_run = Vulcan::Run.first(id: @params[:run_id], workspace_id: @params[:workspace_id])
    if workflow_run
      response = {}
      begin
        config = @remote_manager.read_json_file(workflow_run.run_config_path)
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
