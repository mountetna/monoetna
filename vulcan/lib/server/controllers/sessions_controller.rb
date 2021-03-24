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
    if can_hijack?
      stream_response
    else
      success_json(run_orchestration)
    end
  rescue => e
    Vulcan.instance.logger.log_error(e)
    raise Etna::BadRequest.new(e.message)
  end

  def stream_response
    # We stream the headers back first because
    #   FireFox has a very short timeout for "network-change"
    #   events. If a workflow run takes more than 5s
    #   to return the headers (?), FF will throw a
    #   NetworkError and abort the request.
    # We circumvent this by returning headers immediately and
    #   then stream back the data. The issue will be that
    #   any errors will not be identifiable via status codes.
    @request.env['rack.hijack'].call
    stream = @request.env['rack.hijack_io']

    send_headers(stream)

    begin
      # For FF we have to
      #   continuously send data back, because there appears to be
      #   an additional 10s default timeout such that
      #   if no data is returned within a 10s window, a request
      #   is automatically terminated.
      # Sending a null byte back allows the JSON parser in the client
      #   to still parse the payload without modification.
      # The \n flushes the stream and sends the null byte to the client.
      stream.write("\n")

      Thread.new do
        while !stream.closed?
          stream.write("\n")
          sleep(5)
        end
      end

      stream.write(run_orchestration.to_json)
    rescue => e
      Vulcan.instance.logger.log_error(e)
      # Since we've already sent the status code back as 200,
      #   we'll have to identify errors in the client by
      #   using the error key in the response JSON.
      stream.write({
        error: e.message,
        type: e.class.name
      }.to_json)
    ensure
      stream.close
      @response.close
    end
  end

  def run_orchestration
    orchestration = session.orchestration
    if orchestration.nil?
      raise Etna::NotFound.new("Workflow by the name #{session.workflow_name} could not be found.")
    end

    begin
      orchestration.run_until_done!(storage, token)
    rescue Vulcan::Orchestration::RunErrors => e
      run_errors = e
    end

    {
      session: session.as_json,
      status: orchestration_status(orchestration, run_errors),
      outputs: step_status_json(
        step_name: :primary_outputs,
        bt: orchestration.build_target_for(:primary_outputs)
      ).slice(:status, :downloads)
    }
  end

  def orchestration_status(orchestration, run_errors)
    build_target_cache = {}
    orchestration.unique_paths.map do |paths|
      paths.map do |step_name|
        step = orchestration.workflow.find_step(step_name)
        next nil if step.nil?
        bt = orchestration.build_target_for(step_name, build_target_cache)
        step_status_json(
          step_name: step.id,
          bt: bt,
          run_errors: run_errors,
          ui_output: step.ui_output_name)
      end.select { |v| v }
    end
  end

  def step_status_json(step_name:, bt:, run_errors: nil, ui_output: false)
    {
      name: step_name,
      status: run_errors&.include?(bt) ? 'error' : step_status(bt, ui_output),
      message: run_errors&.include?(bt) ? run_errors.message_for_build_target(bt) : nil,
      downloads: step_has_downloads?(bt, ui_output) ? bt.build_outputs.map do |output_name, sf|
        [
            output_name,
            storage.data_url(project_name: sf.project_name, cell_hash: sf.cell_hash, data_filename: sf.data_filename),
        ]
      end.to_h : nil
    }
  end

  def step_status(bt, ui_output)
    # UI Sinks (plot, download data, etc), will not have outputs
    #   defined, only inputs. So we calculate their
    #   status based on their input file availability.
    return bt.is_buildable?(storage) ? 'complete' : 'pending' if ui_output

    bt.is_built?(storage) ? 'complete' : 'pending'
  end

  def step_has_downloads?(bt, ui_output)
    # UI Sinks (plot, download data, etc), will not have outputs
    #   defined, only inputs. So there are no downloads.
    return false if ui_output

    bt.is_built?(storage)
  end

  def send_headers(stream)
    headers = [
      "HTTP/1.1 200 OK",
      "Content-Type: application/json"
    ]
    stream.write(headers.map { |header| header + "\r\n" }.join)
    stream.write("\r\n")
    stream.flush
  rescue
    stream.close
    raise
  end

  def can_hijack?
    @request.env['rack.hijack?']
  end

  def token
    # For development, we can inject a production token
    #   via config.yml, to talk directly to production services.
    #   This should reduce duplication of data.
    # Note: You'll need to configure the Magma host
    #   in config.yml :development to also point to
    #   production Magma.
    return Vulcan.instance.config(:archimedes_token) || @user.token if :development == Vulcan.instance.environment

    @user.token
  end
end

