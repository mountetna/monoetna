require 'json'
require_relative './vulcan_controller'

class SessionsController < Vulcan::Controller
  def session
    begin
      @the_session ||= Session.from_json(@params.slice(:key, :project_name, :workflow_name, :inputs))
    rescue => e
      raise Etna::BadRequest.new(e.message)
    end
  end

  def orchestration
    session.orchestration
  end

  def scheduler
    orchestration.scheduler
  end

  def workflow
    orchestration.workflow
  end

  def create
    begin
      @the_session = Session.from_json(@params.slice(:project_name, :workflow_name))
    rescue => e
      raise Etna::BadRequest.new(e.message)
    end

    janus_client = Etna::Clients::Janus.new(token: @user.token, host: Vulcan.instance.config(:janus)[:host])
    token = janus_client.generate_token('task', project_name: @params[:project_name], read_only: true)

    success_json({
      session: session.as_json,
      token: token,
    })
  end

  def submit
    orchestration.load_json_inputs!(storage)
    scheduler.schedule_more!(orchestration: orchestration, token: token, storage: storage)
    success_json(status_payload)
  rescue => e
    Vulcan.instance.logger.log_error(e)
    raise Etna::BadRequest.new(e.message)
  end

  def status
    success_json(status_payload)
  end

  def status_payload(build_target_cache: {})
    {
      session: session.as_json,
      status: all_steps_status(build_target_cache: build_target_cache),
      outputs: step_status_json(step_name: :primary_outputs, build_target_cache: build_target_cache).slice(:status, :downloads)
    }
  end

  def all_steps_status(build_target_cache:)
    orchestration.serialized_step_path.map do |paths|
      paths.map do |step_name|
        step = orchestration.workflow.find_step(step_name)
        next nil if step.nil?
        step_status_json(step_name: step.id, build_target_cache: build_target_cache)
      end.select { |v| v }
    end
  end

  def step_status_json(step_name:, build_target_cache:)
    step = workflow.find_step(step_name)
    ui_output = step&.ui_output_name
    bt = orchestration.build_target_for(step_name, build_target_cache)

    scheduler.status(storage: storage, build_target: bt, step: step).update({
      name: step_name,
      downloads: step_has_downloads?(bt, ui_output) ? bt.build_outputs.map do |output_name, sf|
        [
          output_name,
          storage.data_url(project_name: sf.project_name, cell_hash: sf.cell_hash, data_filename: sf.data_filename),
        ]
      end.to_h : nil
    })
  end

  def step_has_downloads?(bt, ui_output)
    # UI Sinks (plot, download data, etc), will not have outputs
    #   defined, only inputs. So there are no downloads.
    return false if ui_output
    bt.is_built?(storage)
  end

  def token
    # For development, we can inject a production token
    #   via config.yml, to talk directly to production services.
    #   This should reduce duplication of data.
    # Note: You'll need to configure the Magma host
    #   in config.yml :development to also point to
    #   production Magma.
    return Vulcan.instance.config(:archimedes_token) || @user.token if :development == Vulcan.instance.environment

    raise Etna::Unauthorized.new("Executing and polling data requires a task token.") unless @user.task?

    @user.token
  end
end

