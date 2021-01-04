require 'open3'

require_relative 'controller'
require_relative '../../etls/redcap/redcap_etl_script_runner'

require 'stringio'

# Add this Kernel method to capture
#   the stdout messages from the ETL
#   process. We'll stream those back
#   to the consumer, since loading
#   can be a long process.
# This will probably not be necessary
#   once we migrate to an async job
#   scheduler.
module Kernel

  def capture_stdout
    out = StringIO.new
    $stdout = out
    yield
    return out
  ensure
    $stdout = STDOUT
  end

end

class JobController < Polyphemus::Controller
  def submit
    require_params(:project_name, :model_names, :redcap_tokens)

    raise Etna::BadRequest, "redcap_tokens must be an array of tokens." unless @params[:redcap_tokens].is_a?(Array)
    raise Etna::BadRequest, "model_names must be [\"all\"] or an array of model names." unless @params[:model_names].is_a?(Array)
    require 'pry'
    binding.pry
    out = capture_stdout do
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

      records = redcap_etl.run(magma_client: magma_client, commit: true)

      return_data = {
          records: records
      }
    end

    return [200, {}, out.string.split("\n")]
    # yield out
    # return success(return_data.to_json, 'application/json')
  end
end
