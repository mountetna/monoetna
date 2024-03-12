require 'etna'
require 'digest'
require_relative "./vulcan_controller"
require_relative './../../ssh'


# TODO: sort out disk location
WORKFLOW_BASE_DIR = "/app/workflows"
WORKSPACE_BASE_DIR = "/app/workspace"
METIS_DIR_MIRROR = "/app/workspace/{HASH}/metis/"
ALLOWED_DIRECTORIES = ["/app/workspace/"]

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

  def create_workflow
    workflow = Vulcan::WorkflowV2.first(project: @params[:project_name], workflow_name: @params[:workflow_name])
    response = {}
    if workflow
      msg = "Workflow: #{@params[:workflow_name]} for project: #{@params[:project_name]} already exists."
      response = {'Warning': msg}
    else
      begin
        # Make project directory if it doesnt exist
        mkdir_output = @ssh.mkdir(project_dir(@escaped_params[:project_name]))

        # Check if there is a repo in that directory
        repo_local_path = repo_local_path(@escaped_params[:project_name], File.basename(@escaped_params[:repo]))
        if @ssh.dir_exists?("#{repo_local_path}")
          Vulcan.instance.logger.info("Repository already exists at: #{repo_local_path} , skipping clone...")
        else
          @ssh.clone(@escaped_params[:repo], @escaped_params[:branch], repo_local_path)
        end

        # Create obj
        obj = Vulcan::WorkflowV2.create(
          project: @escaped_params[:project_name],
          workflow_name: @escaped_params[:workflow_name],
          author: @escaped_params[:author],
          repo_remote_url: @escaped_params[:repo],
          repo_local_path: repo_local_path,
          created_at: Time.now,
          updated_at: Time.now
        )
        response = {'workflow_id': obj.id, 'workflow_name': obj.workflow_name}
      rescue => e
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
        @ssh.mkdir(workspace_dir)
        @ssh.clone(workflow.repo_local_path, @escaped_params[:branch], workspace_dir)
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

  def workflow_run
    success_json({'it works!': true})
  end



end
