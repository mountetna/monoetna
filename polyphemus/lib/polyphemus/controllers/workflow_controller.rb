require_relative 'controller'
require_relative '../../data_eng/workflow_manifests/manifest'
require_relative '../../data_eng/argo_workflow_manager'

class WorkflowController < Polyphemus::Controller

  def create
    require_params(:workflow_name, :workflow_type)

    unless Polyphemus::WorkflowManifest.from_workflow_name(@params[:workflow_type])
      raise Etna::BadRequest, "There is no such workflow type #{@params[:workflow_type]}"
    end

    config = Polyphemus::Config.create(
      project_name: @params[:project_name],
      workflow_name: @params[:workflow_name],
      workflow_type: @params[:workflow_type],
      config_id: Polyphemus::Config.next_id,
      version_number: 1,
      config: {},
      secrets: {},
      created_at: Time.now,
      updated_at: Time.now
    )

    Polyphemus::RuntimeConfig.create(
      config_id: config.config_id,
      run_interval: 0,
      config: {},
      created_at: Time.now,
      updated_at: Time.now
    )

    event_log(
      event: 'create_workflow',
      message: "created #{@params[:workflow_type]} workflow #{@params[:workflow_name]}"
    )

    success_json(config.as_json)
  end

  def update
    config = Polyphemus::Config.current.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).first

    if !config
      raise Etna::NotFound, "No such config #{@params[:config_id]} configured for project #{@params[:project_name]}"
    end
    update = @params.slice(*(config.columns - [:id, :project_name, :config_id]))
    manifest = Polyphemus::WorkflowManifest.from_workflow_name(config.workflow_type)

    if update[:secrets]
      error = manifest.validate_secrets(update[:secrets])
      raise Etna::BadRequest, error if error
      update[:secrets] = config.secrets.merge(update[:secrets])
    end

    if update[:config]
      unless (errors = manifest.validate_config(update[:config])).empty?
        raise Etna::BadRequest, "Invalid configuration for workflow \"#{config.workflow_type}\"\n#{
          errors.map do |error|
            JSONSchemer::Errors.pretty(error)
          end.join("\n")
        }"
      end
      # Here we created a brand new config object, since we have new config values
      update[:version_number] = config.version_number + 1
      update[:updated_at] = Time.now
      update[:created_at] = Time.now
      new_config = Polyphemus::Config.create(config.as_json.merge(update).merge(secrets: config.secrets))
      return success_json(new_config)
    else
      # Here we just update the existing config object
      update[:updated_at] = Time.now
      config.update(update)
    end

    event_log(
      event: 'update_workflow',
      message: "updated workflow #{config.workflow_name} to version #{update[:version_number]}"
    )

    success_json(config.as_json)
  end

  def get_workflows
    return success_json(Polyphemus::WorkflowManifest.list.map(&:as_json))
  end

  def list_by_id
    if @params[:version]
      config = Polyphemus::Config.where(
        project_name: @params[:project_name],
        config_id: @params[:config_id],
        version_number: @params[:version]
      ).first
    else
      config = Polyphemus::Config.current.where(
        project_name: @params[:project_name],
        config_id: @params[:config_id]
      ).first
    end
    success_json(config.as_json)
  end

  def list_all
    if @params[:workflow_type] && !Polyphemus::WorkflowManifest.from_workflow_name(@params[:workflow_type])
      raise Etna::BadRequest, "There is no such workflow type #{@params[:workflow_type]}"
    end

    config_query = Polyphemus::Config.current

    if @params[:workflow_type]
      config_query = config_query.where(workflow_type: @params[:workflow_type])
    end

    configs = config_query.all

    success_json(configs: configs.map(&:as_json))
  end

  def update_run
    require_params(:run_id)
    run = Polyphemus::Run.where(run_id: @params[:run_id]).first

    update_columns = {
      state: @params[:state],
      orchestrator_metadata: @params[:orchestrator_metadata], 
      updated_at: Time.now
    }.compact

    if @params[:output]
      update_columns[:output] = @params[:append_output] ? 
          [run.output, @params[:output]].compact.join : 
          @params[:output]
    end

    if run
      if run.state
        update_columns[:state] = run.state.merge(update_columns[:state] || {})
      end
      run.update(update_columns)
    else
      # For new runs we still need config_id and version_number
      require_params(:config_id, :version_number)
      update_columns.merge!({
        run_id: @params[:run_id],
        name: @params[:name],
        config_id: @params[:config_id],
        version_number: @params[:version_number],
        created_at: Time.now
      })
      run = Polyphemus::Run.create(update_columns)
    end
    success_json(run.as_json)
  end

  def get_run
    run = Polyphemus::Run.where(
      run_id: @params[:run_id]
    ).first

    raise Etna::NotFound, "No such run #{@params[:run_id]}" unless run
    success_json(run.as_json)
  end

  def run_output
    run = Polyphemus::Run.where(
      run_id: @params[:run_id]
    ).first

    raise Etna::NotFound, "No such run #{@params[:run_id]}" unless run
    success_json(output: run.output)
  end

  def runs
    runs = Polyphemus::Run.where(
      config_id: @params[:config_id]
    ).all

    raise Etna::NotFound, "No such config #{@params[:config_id]}" if runs.empty?
    success_json(runs: runs.map(&:summary_hash))
  end

  def get_previous_state
    require_params(:config_id, :state)

    runs = Polyphemus::Run.where(
      config_id: @params[:config_id],
    ).exclude(state: nil).order(Sequel.desc(:created_at)).all

    return success_json({}) if runs.empty?

    filtered_state = @params[:state].to_h do |state_key|
      result = @params[:collect] ?
        runs.map do |run|
          run.state.dig(state_key)
        end.compact : runs.find do |run|
          run.state.has_key?(state_key)
        end&.state&.dig(state_key)
      [ state_key, result ]
    end
    success_json(filtered_state)
  end

  def update_runtime_config
    require_params(:config_id)

    runtime_config = Polyphemus::RuntimeConfig.for_config(@params[:config_id])

    raise Etna::NotFound, "No runtime config found for config_id #{@params[:config_id]}" unless runtime_config

    config = runtime_config.workflow_config

    raise Etna::NotFound, "Cannot find a config for project #{@params[:project_name]} with config_id #{@params[:config_id]}" unless config 

    if @params[:config]
      manifest = Polyphemus::WorkflowManifest.from_workflow_name(config.workflow_type)
      errors = manifest.validate_runtime_config(@params[:config])
      raise Etna::BadRequest, errors.join(', ') if errors.any?
    end
    
    update_columns = {
      config_id: @params[:config_id],
      run_interval: @params[:run_interval],
      config: @params[:config],
      disabled: @params[:disabled],
      updated_at: Time.now
    }.compact

    runtime_config.update(update_columns)
    success_json(runtime_config.as_json)
  end

  def get_runtime_config
    runtime_config = Polyphemus::RuntimeConfig.for_config(@params[:config_id])

    return success_json(runtime_config.as_json) if runtime_config

    raise Etna::NotFound, "No such config #{@params[:config_id]}"
  end

  def run_once
    require_param(:config)
    config = Polyphemus::Config.current.where(
          project_name: @params[:project_name],
          config_id: @params[:config_id]
    ).first

    raise Etna::NotFound, "Cannot find a config for project #{@params[:project_name]} with config_id #{@params[:config_id]}" unless config 

    unless Polyphemus::WorkflowManifest.from_workflow_name(config.workflow_type)
      raise Etna::BadRequest, "There is no such workflow type #{config.workflow_type}"
    end

    manifest = Polyphemus::WorkflowManifest.from_workflow_name(config.workflow_type)
    errors = manifest.validate_runtime_config(@params[:config])
    raise Etna::BadRequest, errors.join(', ') if errors.any?

    begin
      run = Polyphemus::ArgoWorkflowManager.submit_workflow(config)
    rescue StandardError => e
      raise Etna::BadRequest, "Failed to submit Argo workflow: #{e.message}"
    end

    event_log(
      event: 'launch_workflow',
      message: "ran workflow #{config.workflow_name}"
    )
    success_json(run.as_json)
  end

  def revisions
    etl_configs = Polyphemus::Config.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).reverse_order(:version_number).all

    success_json(etl_configs.map(&:to_revision))
  end

  def status
    # Get all unique workflows created by the user (unique workflow_names)
    configs = Polyphemus::Config.current.where(project_name: @params[:project_name]).select(:config_id, :workflow_type, :workflow_name).distinct(:workflow_name).all

    configs = configs.reject do |config|
      config.runtime_config.disabled
    end

    statuses = configs.map do |config|
      # Get the latest run for this config
      # All workflows that have been launched should have a run object
      run, prev_run = Polyphemus::Run.where(config_id: config.config_id).order(Sequel.desc(:created_at)).limit(2).all

      status = {
          workflow_type: config.workflow_type,
          workflow_name: config.workflow_name,
          config_id: config.config_id,
          pipeline_state: nil,
          pipeline_finished_at: nil
      }

      next status.merge(pipeline_state: "pending") if config.runtime_config&.should_run?

      next status unless run

      if run.is_finished?
        next status.merge(
          pipeline_state: run.status,
          pipeline_finished_at: run.finished_at
        )
      else
        next status.merge(
          pipeline_state: "running",
          pipeline_finished_at: prev_run&.is_finished? ? prev_run.finished_at : nil,
        )
      end
    end
    success_json(statuses)
  end
end
