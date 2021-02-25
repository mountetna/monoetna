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

    orchestration.run_until_done!(storage)

    build_target_cache = {}
    success_json({
        session: session.as_json,
        status: orchestration.unique_paths.map do |paths|
          paths.map do |step_name|
            step = orchestration.workflow.find_step(step_name)
            next nil if step.nil?
            step_status_json(step.id, orchestration, build_target_cache)
          end.select { |v| v }
        end,
        outputs: step_status_json(
          :primary_outputs,
          orchestration).slice(:status, :downloads),
    })
  rescue => e
    raise Etna::BadRequest.new(e.message)
  end

  def step_status_json(step_name, orchestration, build_target_cache = {})
    bt = orchestration.build_target_for(step_name, build_target_cache)

    {
        name: step_name,
        status: orchestration.build_target_status(build_target: bt, storage: storage),
        message: orchestration.build_target_has_error?(bt) ? orchestration.build_target_error(bt) : nil,
        downloads: bt.is_built?(storage) ? bt.build_outputs.map do |output_name, sf|
          [
              output_name,
              storage.data_url(project_name: sf.project_name, cell_hash: sf.cell_hash, data_filename: sf.data_filename),
          ]
        end.to_h : nil
    }
  end
end

