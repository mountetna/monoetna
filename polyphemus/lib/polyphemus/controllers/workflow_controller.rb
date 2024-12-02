require_relative 'controller'

class WorkflowController < Polyphemus::Controller
  def get_config
    config_id = @params[:config_id]
    project_name = @params[:project_name]

    config = Polyphemus::Config.where(project_name: project_name, config_id: config_id).first

    unless config
      raise Etna::NotFound, "No such config #{config_id} in project #{project_name}"
    end

    success_json(config.as_json)
  end

  def get_workflow_state
    run_id = @params[:run_id]
    project_name = @params[:project_name]

    run = Polyphemus::Run.where(run_id: run_id).first

    unless run && run.config.project_name == project_name
      raise Etna::NotFound, "No such run #{run_id} in project #{project_name}"
    end

    success_json(run.state)
  end

  def update_workflow_state
    require_params(:state)

    run_id = @params[:run_id]
    project_name = @params[:project_name]
    new_state = @params[:state]

    run = Polyphemus::Run.where(run_id: run_id).first

    unless run && run.config.project_name == project_name
      raise Etna::NotFound, "No such run #{run_id} in project #{project_name}"
    end

    run.update(state: new_state, updated_at: Time.now)

    success_json(run.as_json)
  end

  def write_run_metadata
    require_params(:run_id, :meta_data, :comment)

    run_id = @params[:run_id]
    project_name = @params[:project_name]
    meta_data = @params[:meta_data]
    comment = @params[:comment]

    run = Polyphemus::Run.where(run_id: run_id).first

    unless run && run.config.project_name == project_name
      raise Etna::NotFound, "No such run #{run_id} in project #{project_name}"
    end

    Polyphemus::RunMetadata.create(
      run_id: run_id,
      config_id: run.config_id,
      meta_data: meta_data,
      comment: comment,
      created_at: Time.now,
      updated_at: Time.now
    )

    success_json({ message: "Metadata written successfully" })
  end
end
