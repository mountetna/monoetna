require_relative 'etl_job'
require_relative '../../etls/redcap/redcap_etl_script_runner'

class RedcapLoaderJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithLogger

  def initialize(token, config, runtime_config)
    super
    @workflow_config_id = config['config_id']
    @workflow_version = config['version_number']
  end

  def project_name
    @config['project_name']
  end

  def secrets
    @config['secrets']
  end

  def model_names
    @runtime_config['config']['model_names']
  end

  def config
    @config['config']
  end

  def mode
    @runtime_config['config']['mode']
  end

  def commit?
    !!@runtime_config['config']['commit']
  end

  def pre(context)
    true
  end

  # Process method containing the main File Discovery ETL logic
  def process(context)
    redcap_etl = Polyphemus::RedcapEtlScriptRunner.new(
      project_name: project_name,
      model_names: model_names,
      redcap_tokens: secrets['redcap_tokens'],
      dateshift_salt: Polyphemus.instance.config(:dateshift_salt),
      redcap_host: Polyphemus.instance.config(:redcap)[:host],
      magma_host: Polyphemus.instance.config(:magma)[:host],
      mode: mode,
      config: config
    )

    redcap_etl.run(magma_client: magma_client, commit: commit?, logger: $stdout)
  end

  # Post-condition method to update the number of files to update in the DB
  def post(context)
    polyphemus_client.update_run(project_name, run_id, {
      name: workflow_name,
      config_id: @workflow_config_id,
      version_number: @workflow_version,
      state: { }
    })
  end
end
