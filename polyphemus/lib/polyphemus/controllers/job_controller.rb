require_relative 'controller'
require_relative '../../etls/redcap/redcap_etl_script_runner'

class JobController < Polyphemus::Controller
  def submit
    require_params(:project_name, :model_names, :redcap_tokens)

    raise Etna::BadRequest, "redcap_tokens must be an array of tokens." unless @params[:redcap_tokens].is_a?(Array)
    raise Etna::BadRequest, "model_names must be [\"all\"] or an array of model names." unless @params[:model_names].is_a?(Array)

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
    return success(return_data.to_json, 'application/json')
  end
end
