require 'etna'
require 'digest'
require_relative "./vulcan_controller"
require_relative "./../../path"
require_relative './../../remote_manager'
require_relative './../../snakemake_remote_manager'
require_relative './../../snakemake_command'


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
         @params[:project_name])
      ).all.map do |w|
        w.to_hash
      end
  end

  def update_workflow
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
        @remote_manager.clone(workflow.repo_remote_url, @escaped_params[:branch], workspace_dir)
        @remote_manager.checkout_version(workspace_dir, @escaped_params[:git_version])
        @remote_manager.mkdir(Vulcan::Path.workspace_tmp_dir(workspace_dir))
        @remote_manager.mkdir(Vulcan::Path.workspace_output_path(workspace_dir))
        @remote_manager.upload_dir(Vulcan::Path::SNAKEMAKE_UTILS_DIR, workspace_dir, true)
        config = @remote_manager.read_yaml_file(Vulcan::Path.snakemake_config(workspace_dir))
        obj = Vulcan::Workspace.create(
          workflow_id: workflow.id,
          name: @params[:workspace_name],
          target_mapping: @snakemake_manager.generate_target_mapping(workspace_dir, config),
          path: workspace_dir,
          user_email: @user.email,
          created_at: Time.now,
          updated_at: Time.now
        )
        response = {
          workspace_id: obj.id,
          workflow_id: obj.workflow_id,
          vulcan_config: @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(workspace_dir)),
          dag: @snakemake_manager.get_dag(workspace_dir)
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
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    last_config = Vulcan::Config.where(workspace_id: workspace.id).order(Sequel.desc(:created_at)).first
    last_job_status = Vulcan::Run.where(workspace_id: workspace.id, config_id: last_config.id).first
    response = {
      workspace_id: workspace.id,
      workflow_id: workspace.workflow_id,
      vulcan_config: @remote_manager.read_yaml_file(Vulcan::Path.vulcan_config(workspace.path)),
      dag: @snakemake_manager.get_dag(workspace_dir),
      last_config: last_config.to_hash,
      last_job_status: Vulcan::Run.where(workspace_id: workspace.id).order(Sequel.desc(:created_at)).first.to_hash
    }
    success_json(response)
  end

  def get_dag
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    begin
      success_json({dag: @snakemake_manager.get_dag(workspace.path)})
    rescue => e
      Vulcan.instance.logger.log_error(e)
      raise Etna::BadRequest.new(e.message)
    end
  end

  def save_config
    workspace = Vulcan::Workspace.first(id: @params[:workspace_id])
    unless workspace
      msg = "Workspace for project: #{@params[:project_name]} does not exist."
      raise Etna::BadRequest.new(msg)
    end
    begin
      request_params = @params[:params]
      params_hash = Digest::MD5.hexdigest(request_params.to_json)
      config = Vulcan::Config.first(workspace_id: workspace.id, hash: params_hash)
      unless config
        config_path = Vulcan::Path.workspace_config_path(workspace.path, params_hash)
        @remote_manager.write_file(config_path, request_params.to_json)
        config = Vulcan::Config.create(
          workspace_id: workspace.id,
          path: config_path,
          hash: params_hash,
          created_at: Time.now,
          updated_at: Time.now
        )
      end
      command = Vulcan::Snakemake::CommandBuilder.new
      command.target_meta = {
        mapping: workspace.target_mapping.keys,
        provided_params: request_params.keys.map(&:to_s),
        available_files: @remote_manager.list_files(Vulcan::Path.workspace_output_path(workspace.path)),
      }
      command.options = {
        config_path: config.path,
        profile_path: Vulcan::Path.profile_dir(workspace.path),
        dry_run: true,
      }
      # TODO: do something here with the output
      #@snakemake_manager.run_snakemake(workspace.path, command.build)
      success_json({config_id: config.id})
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
      command = Vulcan::Snakemake::CommandBuilder.new
      command.target_meta = {
        mapping: workspace.target_mapping.transform_values { |v| v.transform_keys(&:to_sym) }, # TODO: maybe change this
        provided_params: @remote_manager.read_json_file(config.path).keys,
        available_files: @remote_manager.list_files(Vulcan::Path.workspace_output_path(workspace.path)).map { |file| "output/#{file}" },
      }
      command.options = {
        config_path: config.path,
        profile_path: Vulcan::Path.profile_dir(workspace.path),
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
      job_id_hash = @snakemake_manager.parse_log_for_slurm_ids(workflow_run.jobs, workflow_run.log_path)
      # TODO: parse the log file for the jobs or do something else
      response = @snakemake_manager.query_sacct(workflow_run.slurm_run_uuid, job_id_hash)
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

end
