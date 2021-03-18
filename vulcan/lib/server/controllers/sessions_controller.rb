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

  def submit
    orchestration.load_json_inputs!(storage)
    scheduler.schedule_more!(orchestration: orchestration, token: @user.token, storage: storage)
    success_json(status_payload)
    # stream.write(run_orchestration.to_json)
  rescue => e
    Vulcan.instance.logger.log_error(e)
    raise Etna::BadRequest.new(e.message)
  end

  def status_payload(build_target_cache: {})
    {
      session: session.as_json,
      status: all_steps_status(build_target_cache: build_target_cache),
      outputs: step_status_json(step_name: :primary_outputs, build_target_cache: build_target_cache).slice(:status, :downloads)
    }
  end

  def all_steps_status(build_target_cache:)
    orchestration.unique_paths.map do |paths|
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

    scheduler.status(storage: storage, build_target: bt).update({
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
end

