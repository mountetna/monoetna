require_relative 'controller'

class EtlController < Polyphemus::Controller
  def jobs
    return success_json(Polyphemus::Job.list.map(&:as_json))
  end

  def list
    success_json(
      Polyphemus::EtlConfig.exclude(archived: true).where(project_name: @params[:project_name]).all.map(&:as_json)
    )
  end

  def revisions
    etl_configs = Polyphemus::EtlConfig.where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).reverse_order(:updated_at).all.sort_by{|e| e.archived ? 1 : 0 }

    success_json(etl_configs.map(&:to_revision))
  end

  def output
    etl_config = Polyphemus::EtlConfig.exclude(archived: true).where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).first

    raise Etna::FileNotFound, "No such etl #{@params[:name]} configured for project #{@params[:project_name]}" unless etl_config

    success_json(output: etl_config.output)
  end

  def update
    etl_configs = Polyphemus::EtlConfig.exclude(archived: true).where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).all

    if etl_configs.empty?
      raise Etna::FileNotFound, "No such etl #{@params[:name]} configured for project #{@params[:project_name]}"
    end

    if etl_configs.length > 1
      begin
        # repair broken configs
        *bad_etl_configs, etl_config = etl_configs.sort_by(&:updated_at)
        Polyphemus::EtlConfig.where(id: bad_etl_configs.map(&:id)).update(archived: true)
      end
    else
      etl_config = etl_configs.first
    end

    update = @params.slice(*(etl_config.columns - [:id, :project_name, :name]))

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
      raise Etna::BadRequest, "Invalid configuration for etl \"#{etl_config.etl}\"" unless etl_config.validate_config(update[:config])

      new_etl_config = Polyphemus::EtlConfig.create(etl_config.as_json.merge(update).merge(secrets: etl_config.secrets))
      etl_config.modified!(:updated_at)
      etl_config.update(archived: true, run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      return success_json(new_etl_config)
    else
      etl_config.update(update)
    end

    success_json(etl_config.as_json)
  end

  def create
    require_params(:job_type)

    etl_configs = Polyphemus::EtlConfig.exclude(archived: true).where(
      project_name: @params[:project_name],
      name: @params[:name]
    ).all

    unless etl_configs.empty?
      raise Etna::BadRequest, "There is already an etl #{@params[:name]} configured for project #{@params[:project_name]}"
    end

    unless Polyphemus::Job.from_name(@params[:job_type])
      raise Etna::BadRequest, "There is no such job type #{@params[:job_type]}"
    end

    etl_config = Polyphemus::EtlConfig.create(
      project_name: @params[:project_name],
      name: @params[:name],
      etl: @params[:job_type],
      config: {},
      secrets: {},
      params: {},
      run_interval: Polyphemus::EtlConfig::RUN_NEVER
    )

    success_json(etl_config.as_json)
  end
end
