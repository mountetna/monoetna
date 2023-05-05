require_relative 'controller'

class EtlController < Polyphemus::Controller
  def jobs
    return success_json(Polyphemus::Job.list.map(&:as_json))
  end

  def list
    configs = Polyphemus::EtlConfig.current.where(project_name: @params[:project_name]).all
    success_json(configs.map(&:as_json))
  end

  def list_all
    if @params[:job_type] && !Polyphemus::Job.from_name(@params[:job_type])
      raise Etna::BadRequest, "There is no such job type #{@params[:job_type]}"
    end

    config_query = Polyphemus::EtlConfig.current

    if @params[:job_type]
      config_query = config_query.where(etl: @params[:job_type])
    end

    configs = config_query.all

    success_json(configs: configs.map(&:as_json))
  end

  def revisions
    etl_configs = Polyphemus::EtlConfig.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).reverse_order(:version_number).all

    success_json(etl_configs.map(&:to_revision))
  end

  def output
    etl_config = Polyphemus::EtlConfig.current.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).first

    raise Etna::NotFound, "No such etl #{@params[:config_id]} configured for project #{@params[:project_name]}" unless etl_config

    success_json(output: etl_config.output)
  end

  def update
    etl_config = Polyphemus::EtlConfig.current.where(
      project_name: @params[:project_name],
      config_id: @params[:config_id]
    ).first

    if !etl_config
      raise Etna::NotFound, "No such etl #{@params[:config_id]} configured for project #{@params[:project_name]}"
    end

    update = @params.slice(*(etl_config.columns - [:id, :project_name, :config_id]))

    if update[:secrets]
      error = etl_config.validate_secrets(update[:secrets])

      raise Etna::BadRequest, error if error

      update[:secrets] = etl_config.secrets.merge(update[:secrets])
    end

    if update[:params]
      errors = etl_config.validate_params(update[:params])

      unless errors.empty?
        raise Etna::BadRequest, errors.join('; ').capitalize
      end
      update[:params] = etl_config.params.merge(update[:params])
    end

    if update[:config]
      unless (errors = etl_config.validate_config(update[:config])).empty?
        raise Etna::BadRequest, "Invalid configuration for etl \"#{etl_config.etl}\"\n#{
          errors.map do |error|
            JSONSchemer::Errors.pretty(error)
          end.join("\n")
        }"
      end

      update[:version_number] = etl_config.version_number + 1

      new_etl_config = Polyphemus::EtlConfig.create(etl_config.as_json.merge(update).merge(secrets: etl_config.secrets))
      return success_json(new_etl_config)
    else
      etl_config.update(update)
    end

    success_json(etl_config.as_json)
  end

  def create
    require_params(:name, :job_type)

    etl_configs = Polyphemus::EtlConfig.current.where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).all

    unless Polyphemus::Job.from_name(@params[:job_type])
      raise Etna::BadRequest, "There is no such job type #{@params[:job_type]}"
    end

    etl_config = Polyphemus::EtlConfig.create(
      project_name: @params[:project_name],
      name: @params[:name],
      etl: @params[:job_type],
      config_id: Polyphemus::EtlConfig.next_id,
      version_number: 1,
      config: {},
      secrets: {},
      params: {},
      run_interval: Polyphemus::EtlConfig::RUN_NEVER
    )

    success_json(etl_config.as_json)
  end
end
