require 'action_controller'

require_relative 'controller'
require_relative '../../etls/redcap/redcap_etl_script_runner'


class JobController < Polyphemus::StreamingController
  def submit
    require_params(:project_name, :model_names, :redcap_tokens)

    raise Etna::BadRequest, "redcap_tokens must be an array of tokens." unless @params[:redcap_tokens].is_a?(Array)
    raise Etna::BadRequest, "model_names must be \"all\" or an array of model names." unless @params[:model_names].is_a?(Array) || "all" == @params[:model_names]

    if can_hijack?
      run_loading_in_thread
    else
      # Mostly for testing and as a basic fallback...
      return run_loading_inline
    end
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end

  private

  def run_loading_in_thread
    # Rack hijacking technique found here:
    #   https://blog.chumakoff.com/en/posts/rails_sse_rack_hijacking_api
    @request.env['rack.hijack'].call
    stream = @request.env['rack.hijack_io']

    send_headers(stream)

    Thread.new do
      launch_redcap_process(stream)
    end

    @response.close
  end

  def run_loading_inline
    records = run_redcap_loader(STDOUT)
    success({record: records}.to_json, 'application/json')
  end

  def run_redcap_loader(logger)
    redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
      project_name: @params[:project_name],
      model_names: @params[:model_names],
      redcap_tokens: @params[:redcap_tokens],
      dateshift_salt: Polyphemus.instance.config(:dateshift_salt).to_s,
      redcap_host: Polyphemus.instance.config(:redcap)[:host],
      magma_host: Polyphemus.instance.config(:magma)[:host]
    )

    magma_client = Etna::Clients::Magma.new(
      token: @user.token,
      host: Polyphemus.instance.config(:magma)[:host])

    # TODO: Need to change commit to true before merge + deploy
    redcap_etl.run(magma_client: magma_client, commit: false, logger: logger)
  end

  def launch_redcap_process(stream)
    sse = ActionController::Live::SSE.new(stream, retry: 300, event: "REDCapLoadingProgress")

    run_redcap_loader(sse)
  rescue => e
    sse.write("#{e.message}\n#{e.backtrace}.")
  ensure
    sse.close
  end
end
