require_relative 'controller'
require_relative '../../etls/redcap/redcap_etl_script_runner'

class JobController < Polyphemus::Controller
  def submit
    require_params(:project_name, :model_names, :redcap_token)

    redcap_loader = RedcapEtlScriptRunner.new(
      project_name: @params[:project_name],
      model_names: @params[:model_names],
      redcap_token: @params[:redcap_token]
    )

    magma_client = Etna::Clients::Magma.new(
      token: @user.token,
      host: Polyphemus.instance.config(:magma, environment)[:host])

    records = redcap_loader.run(magma_client: magma_client)

    return_data = {
        records: records
    }
    return success(return_data.to_json, 'application/json')
  end
end
