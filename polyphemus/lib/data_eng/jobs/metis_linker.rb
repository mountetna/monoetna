require_relative 'etl_job'
require_relative '../../etls/metis/loader'

class MetisLinkerJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithLogger

  def project_name
    config['project_name']
  end

  def workflow_config_id
    config['config_id']
  end

  def workflow_version
    config['version_number']
  end

  def bucket_name
    config['config']['bucket_name']
  end

  def pre(context)
    context[:start_time] = Time.at(fetch_last_scan).to_datetime.iso8601
    context[:end_time] = Time.now.to_datetime.iso8601
    true
  end

  # Process method containing the main File Discovery ETL logic
  def process(context)
    last_scan = fetch_last_scan

    rules = gnomon_client.project_rules(project_name).rules

    project_def = magma_client.retrieve(project_name: project_name)

    tail = metis_client.tail_bucket(Etna::Clients::Metis::TailBucketRequest.new(
        project_name: project_name,
        bucket_name: bucket_name,
        type: 'files',
        batch_start: context[:start_time],
        batch_end: context[:end_time]
    ))

    params = runtime_config['config'] || {}

    loader = Metis::Loader.new(config, rules, params, project_def.models)

    update = loader.update_for(tail, metis_client)

    response = magma_client.update_json(update, page_size=100)

    summary = <<EOT
===============================
Upload Summary : #{context[:start_time]} -> #{context[:end_time]} 
Models: #{response.models.model_keys.join(', ')}
Committed to Magma: #{!loader.config.dry_run?}
Autolinked Parent Identifiers: #{loader.config.autolink?}
EOT

    response.models.each do |model_name, model|
      summary += "#{model_name} records updated: #{model.documents.document_keys.join(', ')}\n"
    end

    summary += "==============================="

    puts summary
  end

  # Post-condition method to update the number of files to update in the DB
  def post(context)
    polyphemus_client.update_run(project_name, run_id, {
      name: workflow_name,
      config_id: workflow_config_id,
      version_number: workflow_version,
      state: {
        start_time: context[:start_time],
        end_time: context[:end_time]
      }
    })
  end

  private

  def collect_tails
  end

  # Fetch the last_scan timestamp from the pipeline state using Polyphemus client
  def fetch_last_scan
    begin
      response = polyphemus_client.get_previous_state(
       project_name,
       workflow_config_id,
       state: [:end_time]
      )
      return response['end_time'].to_i
    rescue Etna::Error => e
      return 0
    end
  end
end
