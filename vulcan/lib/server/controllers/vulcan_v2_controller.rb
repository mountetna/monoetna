require 'etna'
require_relative "./vulcan_controller"


# TODO: make this dynamic
WORKFLOW_DIR = "/app/workflows"

class VulcanV2Controller < Vulcan::Controller

  def invoke_ssh_command(command)
    stdout_data = ""
    stderr_data = ""
    exit_status = nil

    Vulcan.instance.ssh.open_channel do |channel|
      channel.exec(command) do |ch, success|
        unless success
          raise "Command execution failed: #{command}"
        end

        channel.on_data do |ch, data|
          stdout_data += data
        end

        channel.on_extended_data do |ch, type, data|
          stderr_data += data
        end

        channel.on_request("exit-status") do |ch, data|
          exit_status = data.read_long
        end
      end
    end

    Vulcan.instance.ssh.loop

    if exit_status != 0
      raise "Command exited with status #{exit_status}: #{command}"
    end

    { stdout: stdout_data, stderr: stderr_data, exit_status: exit_status }
  end

  def clone(repo_url, branch, target_directory)
    # This will be used to clone a repo that will be
    # instantiated as a workflow and then also to clone a repo
    # into a workspace
    #
    # For some control over what can be cloned we can just have a hardcoded list with
    # available repos to be cloned
    command = "git clone -b #{branch} #{repo_url} #{target_directory} "
    invoke_ssh_command(command)
  end

  def clone_workflow
    # This should really only be invoked once
    # Make project directory if it doesnt exist
    require 'pry'; binding.pry
    command = "mkdir -p #{WORKFLOW_DIR}/#{@params[:project_name]}"
    output = invoke_ssh_command(command)
    puts output

    # Clone workflow
    target_dir = "#{WORKFLOW_DIR}/#{@params[:project_name]}"
    out = clone(@params[:repo], @params[:branch], target_dir )
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

  def params
    success_json({'it works!': true})
  end


  def list_workspaces
    success_json({'it works!': true})
  end


end
