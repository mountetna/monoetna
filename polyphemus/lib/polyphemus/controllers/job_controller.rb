require 'open3'
require 'action_controller'

require_relative 'controller'
require_relative '../../etls/redcap/redcap_etl_script_runner'


# Wrap StringIO in case we need to flush old
#   responses as new stdout messages come in?
# Also, we need to wrap StringIO so that
#   #each returns an actual string for the response.
class StringBuffer
  def initialize
    @strs = []
  end

  def <<(value)
    @strs << value
  end

  def write(value)
    @strs << value
  end

  def each_line
    yield @strs.slice!(0)
  end

  def each
    yield @strs.slice!(0)
  end
end

# Add this Kernel method to capture
#   the stdout messages from the ETL
#   process. We'll stream those back
#   to the consumer, since loading
#   can be a long process.
# This will probably not be necessary
#   once we migrate to an async job
#   scheduler.
# module Kernel

#   def capture_stdout
#     out = StringBuffer.new
#     $stdout = out
#     yield
#     return out
#   ensure
#     $stdout = STDOUT
#   end

# end

class JobController < Polyphemus::StreamingController
  def submit
    puts "in submit"
    require_params(:project_name, :model_names, :redcap_tokens)

    raise Etna::BadRequest, "redcap_tokens must be an array of tokens." unless @params[:redcap_tokens].is_a?(Array)
    raise Etna::BadRequest, "model_names must be \"all\" or an array of model names." unless @params[:model_names].is_a?(Array) || "all" == @params[:model_names]

    # response.headers['Content-Type'] = 'text/event-stream'
    # response.headers['Last-Modified'] = Time.now.httpdate
    # require 'pry'
    # binding.pry
    # out = capture_stdout do
    #   puts "standing up the runner"
    #   redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
    #     project_name: @params[:project_name],
    #     model_names: @params[:model_names],
    #     redcap_tokens: @params[:redcap_tokens],
    #     dateshift_salt: Polyphemus.instance.config(:dateshift_salt).to_s,
    #     redcap_host: Polyphemus.instance.config(:redcap)[:host],
    #     magma_host: Polyphemus.instance.config(:magma)[:host]
    #   )

    #   puts "standing up a magma client"
    #   magma_client = Etna::Clients::Magma.new(
    #     token: @user.token,
    #     host: Polyphemus.instance.config(:magma)[:host])

    #   puts "calling run"
    #   redcap_etl.run(magma_client: magma_client, commit: false)
    # end
    # require 'pry'
    # binding.pry


    # @response.headers["Content-Type"] = "text/event-stream"

    # @response.headers["rack.hijack"] = proc do |stream|
    #   Thread.new do
    #     launch_redcap_process(stream)
    #   end
    # end

    # "head" is ActionController specific?
    # head :ok

    @request.env['rack.hijack'].call
    stream = @request.env['rack.hijack_io']

    send_headers(stream)

    Thread.new do
      launch_redcap_process(stream)
    end

    @response.close

    # redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
    #   project_name: @params[:project_name],
    #   model_names: @params[:model_names],
    #   redcap_tokens: @params[:redcap_tokens],
    #   dateshift_salt: Polyphemus.instance.config(:dateshift_salt).to_s,
    #   redcap_host: Polyphemus.instance.config(:redcap)[:host],
    #   magma_host: Polyphemus.instance.config(:magma)[:host]
    # )

    # puts "standing up a magma client"
    # magma_client = Etna::Clients::Magma.new(
    #   token: @user.token,
    #   host: Polyphemus.instance.config(:magma)[:host])

    # logger = StringBuffer.new
    # puts "calling run"
    # redcap_etl.run(magma_client: magma_client, commit: false, logger: STDOUT)
    # # redcap_etl.run(magma_client: magma_client, commit: false, logger: logger) do |logger|
    # #   logger.each_line do |line|
    # #     response.stream.write(line)
    # #   end
    # # end

    # return [200, {'Last-Modified': Time.now.httpdate}, logger]
  rescue => e
    # puts e.message
    puts "#{e.message}\n#{e.backtrace}"
    Polyphemus.instance.logger.log_error(e)
    return failure("#{e.message}\n#{e.backtrace}")
    # command = "./bin/polyphemus run_redcap_loader #{Polyphemus.instance.environment.to_s} #{@params[:project_name]} #{@params[:model_names]} #{@params[:redcap_tokens]} --execute"
    # Open3.popen3(command) {|stdin, stdout, stderr, wait_thr|
    #   pid = wait_thr[:pid]

    #   stdin.close

    #   stdout.each_line do |line|
    #     puts "stdout"
    #     puts pid
    #     puts line
    #     puts wait_thr.value.exitstatus
    #     @response.body << line
    #     if wait_thr.value.exitstatus
    #       Process.kill("KILL", pid)
    #     end
    #   end

    #   stderr.each_line do |line|
    #     puts "stderr"
    #     puts line
    #     puts pid
    #     puts wait_thr.value.exitstatus
    #     @response.body << line
    #     if wait_thr.value.exitstatus
    #       Process.kill("KILL", pid)
    #     end
    #   end
    # }
      # launch_redcap_process do |line, exit_status|
      #   # binding.pry
      #   puts line
      #   puts exit_status
      #   binding.pry if exit_status
      #   return if exit_status
      #   response.stream.write(line)
      # end
    # Polyphemus.instance.logger.info("Job completed.")
    # response.stream.close


    # return [200, {}, SlowStreamer.new]
    # yield out
    # return success(return_data.to_json, 'application/json')
  # ensure
    # response.stream.close
  end

  private

  def launch_redcap_process(stream)
    sse = ActionController::Live::SSE.new(stream, retry: 300, event: "taskProgress")

    redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
      project_name: @params[:project_name],
      model_names: @params[:model_names],
      redcap_tokens: @params[:redcap_tokens],
      dateshift_salt: Polyphemus.instance.config(:dateshift_salt).to_s,
      redcap_host: Polyphemus.instance.config(:redcap)[:host],
      magma_host: Polyphemus.instance.config(:magma)[:host]
    )

    puts "standing up a magma client"
    magma_client = Etna::Clients::Magma.new(
      token: @user.token,
      host: Polyphemus.instance.config(:magma)[:host])

    puts "calling run"
    redcap_etl.run(magma_client: magma_client, commit: false, logger: sse)
  ensure
    sse.close
  end

  def send_headers(stream)
    headers = [
      "HTTP/1.1 200 OK",
      "Content-Type: text/event-stream"
    ]
    stream.write(headers.map { |header| header + "\r\n" }.join)
    stream.write("\r\n")
    stream.flush
  rescue
    stream.close
    raise
  end
end
