require 'etna'
require 'shellwords'
require_relative "./vulcan_controller"


# TODO: sort out disk location
WORKFLOW_DIR = "/app/workflows"
WORKSPACE_DIR = "/app/workspace"
METIS_DIR_MIRROR = "/app/workspace/{HASH}/metis/"

class VulcanV2Controller < Vulcan::Controller

  def create_workflow
    workflow = Vulcan::Workflow[:project => @params[:project_name], :workflow_name => @params[:workflow_name] ]
    if workflow
      puts workflow
    else
      begin
        # Make project directory if it doesnt exist
        command = "mkdir -p #{WORKFLOW_DIR}/#{@escaped_params[:project_name]}"
        mkdir_output = invoke_ssh_command(command)

        # Clone workflow
        target_dir = "#{WORKFLOW_DIR}/#{@escaped_params[:project_name]}"
        command = "git clone -b #{@escaped_params[:branch]} #{@escaped_params[:repo]} #{target_dir}"
        clone_output = invoke_ssh_command(command)

        # Create obj
        Vulcan::Workflow.create(
          project_name: @escaped_params[:project_name],
          workflow_name: @escaped_params[:workflow_name],
          author: @escaped_params[:author],
          repository_url: @escaped_params[:repository_url],
          created_at: Time.now,
          updated_at: Time.now
        )
      rescue => e
        Vulcan.instance.logger.log_error(e.to_s)
        raise Etna::BadRequest.new(e.to_s)
      end
    end
    success_json({'it works!': true})
  end

  def update_workflow
    success_json({'it works!': true})
  end

  def list_workflows
    # require 'pry'; binding.pry
    command = "ls #{WORKFLOW_DIR}/#{config_json[:project_name]}/"
    out = invoke_ssh_command(command)
    puts out
    success_json({'it works!': true})
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
