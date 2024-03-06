require 'etna'
require 'shellwords'
require_relative "./vulcan_controller"
require_relative './../../ssh'


# TODO: sort out disk location
WORKFLOW_DIR = "/app/workflows"
WORKSPACE_DIR = "/app/workspace"
METIS_DIR_MIRROR = "/app/workspace/{HASH}/metis/"

class VulcanV2Controller < Vulcan::Controller

  def initialize(request, action = nil)
    super
    @ssh = Vulcan::SSH.new(Vulcan.instance.ssh)
  end

  def create_workflow
    workflow = Vulcan::WorkflowV2.first(project: @params[:project_name], workflow_name: @params[:workflow_name])
    msg = ''
    if workflow
      msg = "Workflow: #{@params[:workflow_name]} for project: #{@params[:project_name]} already exists."
    else
      begin
        # Make project directory if it doesnt exist
        mkdir_output = @ssh.mkdir("#{WORKFLOW_DIR}/#{@escaped_params[:project_name]}")

        # Check if there is a repo in that directory
        target_dir = "#{WORKFLOW_DIR}/#{@escaped_params[:project_name]}/#{File.basename(@escaped_params[:repo])}"
        unless @ssh.dir_exists?("#{target_dir}")
          @ssh.clone(@escaped_params[:branch], @escaped_params[:repo], target_dir)
        end

        # Create obj
        Vulcan::WorkflowV2.create(
          project: @escaped_params[:project_name],
          workflow_name: @escaped_params[:workflow_name],
          author: @escaped_params[:author],
          repository_url: @escaped_params[:repo],
          created_at: Time.now,
          updated_at: Time.now
        )
        msg = "Workflow: #{@params[:workflow_name]} successfully cloned and created."
      rescue => e
        Vulcan.instance.logger.log_error(e)
        raise Etna::BadRequest.new(e.message)
      end
    end
    success_json({'info': msg})
  end

  def update_workflow
    success_json({'it works!': true})
  end

  def list_workflows
    success_json(
      figures: Vulcan::WorkflowV2.where(
        project: @params[:project_name]
      ).all.map do |f|
        f.to_hash
      end
    )
  end

  def workspace_create
    success_json({'it works!': true})
  end

  def workflow_run
    success_json({'it works!': true})
  end

  def list_workspaces
    success_json({'it works!': true})
  end


end
