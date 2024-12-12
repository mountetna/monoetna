require_relative 'controller'
require_relative '../../data_eng/workflow_manifests/manifest'

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

      update[:version_number] = config.version_number + 1

      new_config = Polyphemus::Config.create(config.as_json.merge(update).merge(secrets: config.secrets))
      return success_json(new_config)
    else
      config.update(update)
    end

    success_json(config.as_json)
  end

  def workflows
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

  def list
    configs = Polyphemus::Config.current.where(project_name: @params[:project_name]).all
    success_json(configs.map(&:as_json))
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

  def revisions
    etl_configs = Polyphemus::Config.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).reverse_order(:version_number).all

    success_json(etl_configs.map(&:to_revision))
  end

  def output
    etl_config = Polyphemus::Config.current.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).first

    raise Etna::NotFound, "No such etl #{@params[:config_id]} configured for project #{@params[:project_name]}" unless etl_config

    success_json(output: etl_config.output)
  end

  def add_output
    etl_config = Polyphemus::Config.current.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).first

    if !etl_config
      raise Etna::NotFound, "No such etl #{@params[:config_id]} configured for project #{@params[:project_name]}"
    end

    output = [
      (@params[:append] ? etl_config.output : nil),
      @params[:output]
    ].compact.join

    etl_config.update(output: output)

    success_json(etl_config.as_json.merge(output: output))
  end

  
  end
