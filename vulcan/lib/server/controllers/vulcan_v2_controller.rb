require 'etna'
require 'digest'
require_relative "./vulcan_controller"
require_relative './../../ssh'


# TODO: sort out disk location
WORKFLOW_BASE_DIR = "/app/vulcan/workflows"
WORKSPACE_BASE_DIR = "/app/vulcan/workspace"
SNAKEMAKE_UTILS_DIR = "/app/snakemake_utils" # This is local
TMPDIR = "/tmp/vulcan"
ALLOWED_DIRECTORIES = [WORKSPACE_BASE_DIR, WORKFLOW_BASE_DIR, TMPDIR]

class VulcanV2Controller < Vulcan::Controller

  def initialize(request, action = nil)
    super
    @ssh = Vulcan::SSH.new(Vulcan.instance.ssh)
  end

  def project_dir(project_name)
    "#{WORKFLOW_BASE_DIR}/#{project_name}"
  end

  def repo_local_path(project_name, repository_name)
    "#{project_dir(project_name)}/#{repository_name}"
  end

  def workspace_dir(project_name, hash)
    "#{WORKSPACE_BASE_DIR}/#{project_name}/#{hash}"
  end

  def workspace_hash(workflow_id, user_email)
    Digest::MD5.hexdigest(workflow_id + Time.now.to_s + user_email)
  end

  def tmp_hash(workflow_name, user_email)
    Digest::MD5.hexdigest(workflow_name + Time.now.to_s + user_email)
  end

  def tmp_dir(tmp_hash)
    "#{TMPDIR}/#{tmp_hash}"
  end

  def metis_mirror_path(workspace_dir)
    "#{workspace_dir}/metis_output/"
  end

  # Admin command complete prior
  def create_repo
    repo_name = File.basename(@escaped_params[:repo_url])
    repo_local_path = repo_local_path(@escaped_params[:project_name], repo_name)
    if @ssh.dir_exists?("#{repo_local_path}")
      msg = "Repo: #{repo_local_path} for project: #{@params[:project_name]} already exists."
      Vulcan.instance.logger.info(msg)
      response = {'Warning': msg}
    else
      begin
        # Make project directory if it doesnt exist
        mkdir_output = @ssh.mkdir(project_dir(@escaped_params[:project_name]))
        # Create repo directory
        @ssh.clone(@escaped_params[:repo_url], @escaped_params[:branch], repo_local_path)
        response = {repo_local_path: repo_local_path, repo_name: repo_name}
      rescue => e
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
    end
    success_json(response)
  end

  # Admin command
  def list_repos
    begin
      project_dir = project_dir(@escaped_params[:project_name])
      success_json({dirs: @ssh.list_dirs(project_dir)})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
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
        tmp_dir = tmp_dir(tmp_hash(@escaped_params[:workflow_name], @user.email))
        @ssh.mkdir(tmp_dir)
        @ssh.clone(@escaped_params[:repo_local_path], @escaped_params[:branch], tmp_dir)
        @ssh.checkout_tag(tmp_dir, @escaped_params[:tag])
        config = @ssh.read_yaml_file("#{tmp_dir}/vulcan_config.yaml")
        # TODO: run a validation on the config
        obj = Vulcan::WorkflowV2.create(
          project: @escaped_params[:project_name],
          workflow_name: @escaped_params[:workflow_name],
          author: @escaped_params[:author],
          repo_remote_url: @ssh.get_repo_remote_url(@escaped_params[:repo_local_path]),
          repo_local_path: @escaped_params[:repo_local_path],
          repo_tag: @escaped_params[:tag],
          config: config.to_json,
          created_at: Time.now,
          updated_at: Time.now
        )
        @ssh.rmdir(tmp_dir, ALLOWED_DIRECTORIES)
        response = {'workflow_id': obj.id, 'workflow_name': obj.workflow_name}
      rescue => e
        @ssh.rmdir(tmp_dir, ALLOWED_DIRECTORIES)
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
      workspace_hash = workspace_hash(workflow.id.to_s, @user.email)
      workspace_dir = workspace_dir(@escaped_params[:project_name], workspace_hash)
      begin
        @ssh.mkdir(workspace_dir)
        # TODO: we probably want to store the name of the "master" branch
        @ssh.clone(workflow.repo_local_path, "main", workspace_dir)
        @ssh.checkout_tag(workspace_dir, workflow.repo_tag)
        #@ssh.mkdir(metis_mirror_path(workspace_dir)) # location of output files that will get mirrored back to metis
        @ssh.upload_dir(SNAKEMAKE_UTILS_DIR, workspace_dir, true)
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
        @ssh.rmdir(workspace_dir, ALLOWED_DIRECTORIES)
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
    workspace = Vulcan::Workspace.first(
      id: @params[:workspace_id],
    )
    begin
      @ssh.run_snakemake(workspace.path, @params[:config])
      # poll for success
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
    # TODO Validation on steps
  end

end
