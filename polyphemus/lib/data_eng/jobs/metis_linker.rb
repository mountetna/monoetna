require_relative 'etl_job'

class MetisLinkerJob < Polyphemus::ETLJob
  include WithEtnaClients

  def initialize(config, runtime_config)
    @config = config
    @runtime_config = runtime_config
    @workflow_config_id = config['config_id']
    @workflow_version = config['version_number']
  end

  def pre
    true
  end

  # Process method containing the main File Discovery ETL logic
  def process
    last_scan = fetch_last_scan

    rules = gnomon_client.rules(config[:project_name])[:rules]

    models = magma_client.retrieve(
      project_name: project_name,
      model_name: 'all',
      record_names: [], attribute_names: "all"
    )

    tail = metis.tail(
        project_name: project_name,
        bucket_name: bucket_name,
        type: 'files',
        batch_start: ran_at,
        batch_end: Time.now
    )

    update = Metis::Loader.new(config, rules, {}, models).update_for(tail, metis_client)

    magma_client.update(update)
  end

  # Post-condition method to update the number of files to update in the DB
  def post
  end

  private

  def project_name
    @config[:project_name]
  end

  def collect_tails
  end

  # Fetch the last_scan timestamp from the pipeline state using Polyphemus client
  def fetch_last_scan
    run = polyphemus_client.get_run(run_id)
    run.state[:last_scan] ? Time.parse(run.state[:last_scan]) : Time.at(0) # Default to epoch if not found
  rescue StandardError => e
    puts "Error fetching pipeline state: #{e.message}"
    Time.at(0)
  end

  # Fetch files from SFTP that are newer than last_scan using regex
  def fetch_files_from_sftp(last_scan)
    files = []
    entries = @sftp_client.list_directory(config[:sftp_directory])

    entries.each do |entry|
      next if entry[:name] == '.' || entry[:name] == '..'

      full_path = File.join(config[:sftp_directory], entry[:name])
      next unless file_matches_criteria?(entry, last_scan)

      md5 = @sftp_client.calculate_md5(full_path)
      files << [full_path, md5] if md5 != 'error'
    end
    files
  rescue StandardError => e
    puts "Error fetching files from SFTP: #{e.message}"
    raise
  end

  # Check if the file matches the required criteria
  def file_matches_criteria?(entry, last_scan)
    matches_type = entry[:longname] =~ @regex
    newer_than_last_scan = Time.at(entry[:attributes].mtime) > last_scan
    matches_type && newer_than_last_scan
  end
end
