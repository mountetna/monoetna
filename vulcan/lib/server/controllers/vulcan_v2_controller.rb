require 'etna'
require 'digest'
require_relative "./vulcan_controller"
require_relative './../../ssh'


# TODO: sort out disk location
WORKFLOW_BASE_DIR = "/app/vulcan/workflows"
WORKSPACE_BASE_DIR = "/app/vulcan/workspace"
SNAKEMAKE_UTILS_DIR = "/app/vulcan/snakemake_utils"
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
    repo_name = File.basename(@escaped_params[:repo])
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
        @ssh.clone(@escaped_params[:repo], @escaped_params[:branch], repo_local_path)
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
  end

  # User
  def publish_workflow
    workflow = Vulcan::WorkflowV2.first(
      project: @params[:project_name],
      workflow_name: @params[:workflow_name],
      tag: @params[:tag]
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
        # Make sure config is valid
        config = @ssh.read_yaml_file("#{tmp_dir}/config.yaml")
        # TODO: parse config
        obj = Vulcan::WorkflowV2.create(
          project: @escaped_params[:project_name],
          workflow_name: @escaped_params[:workflow_name],
          author: @escaped_params[:author],
          repo_remote_url: "some/remote/path", # TODO
          repo_local_path: @escaped_params[:repo_local_path],
          tag: @escaped_params[:tag],
          params: [],
          jobs: [],
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
    workflow = Vulcan::WorkflowV2.first(project: @params[:project_name], workflow_name: @params[:workflow_name])
    response = {}
    if workflow
      workspace_hash = workspace_hash(workflow.id.to_s, @user.email)
      workspace_dir = workspace_dir(@escaped_params[:project_name], workspace_hash)
      begin
        # TODO: finish this up
        @ssh.mkdir(workspace_dir)
        @ssh.clone(workflow.repo_local_path, @escaped_params[:branch], workspace_dir)
        @ssh.mkdir(metis_mirror_path(workspace_dir)) # location of output files that will get mirrored back to metis
        @ssh.upload_dir(SNAKEMAKE_UTILS_DIR, workspace_dir, true)
        require 'pry'; binding.pry
        #TODO: fix this
        obj = Vulcan::Workspace.create(
          workflow_id: workflow.id,
          workspace_dir: workspace_dir,
          repo_branch: @escaped_params[:branch],
          user_email: @user.email,
          hash: workspace_hash,
          created_at: Time.now,
          updated_at: Time.now
        )
        response = {'workspace_id': obj.id, 'workflow_hash': obj.hash_value, workflow_id: obj.workflow_id}
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
      workflows: Vulcan::Workspace.where(
        user_email: @params[:user_email],
        workflow_id: @params[:workflow_id]
      ).all.map do |w|
        w.to_hash
      end
    )
  end

  def retrieve_file_params
    # Use metis client to download files into /hash/inputs/
    # For extremely large
  end

  def workflow_run
  end



end
