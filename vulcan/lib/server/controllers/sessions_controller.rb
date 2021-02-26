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

  def submit
    orchestration = session.orchestration
    if orchestration.nil?
      raise Etna::NotFound.new("Workflow by the name #{session.workflow_name} could not be found.")
    end

    begin
      orchestration.run_until_done!(storage)
    rescue Vulcan::Orchestration::RunErrors => e
      run_errors = e
    end

    success_json({
        session: session.as_json,
        status: orchestration_status(orchestration, run_errors),
        outputs: step_status_json(
          :primary_outputs,
          orchestration.build_target_for(:primary_outputs)).slice(:status, :downloads)
    })
  rescue => e
    Vulcan.instance.logger.log_error(e)
    raise Etna::BadRequest.new(e.message)
  end

  def orchestration_status(orchestration, run_errors)
    build_target_cache = {}
    orchestration.unique_paths.map do |paths|
      paths.map do |step_name|
        step = orchestration.workflow.find_step(step_name)
        next nil if step.nil?
        bt = orchestration.build_target_for(step_name, build_target_cache)
        step_status_json(step.id, bt, run_errors)
      end.select { |v| v }
    end
  end

  def step_status_json(step_name, bt, run_errors=nil)
    {
      name: step_name,
      status: run_errors&.include?(bt) ? 'error' : bt.is_built?(storage) ? 'complete' : 'pending',
      message: run_errors&.include?(bt) ? run_errors.message_for_build_target(bt) : nil,
      downloads: bt.is_built?(storage) ? bt.build_outputs.map do |output_name, sf|
        [
            output_name,
            storage.data_url(project_name: sf.project_name, cell_hash: sf.cell_hash, data_filename: sf.data_filename),
        ]
      end.to_h : nil
    }
  end
end

