require_relative '../../lib/sftp_client'

class FileDiscoveryJob < ETLJob
  include WithEtnaClients

  def initialize(config, secrets)
    super(config, secrets)
    @sftp_client = SftpClient.new(
      host: config[:sftp_host],
      user: config[:sftp_user],
      password: secrets[:sftp_password],
      port: config[:sftp_port] || 22
    )
    @workflow_config_id = config[:config_id]
    @regex = config[:regex]
    @root_dir = config[:root_dir]
  end

  def pre
    true
  end

  # Process method containing the main File Discovery ETL logic
  def process
    last_scan = fetch_last_scan
    files_to_update = fetch_files_from_sftp(last_scan)
    @sftp_client.write_csv(config[:files_to_update_path], files_to_update)
  end

  # Post-condition method to update the number of files to update in the DB
  def post
    file_count = File.readlines(config[:files_to_update_path]).size
    polyphemus_client.update_workflow_state(argo_id, {
      num_files_to_update: file_count,
      updated_at: Time.now,
      modified_at: Time.now,
    })
    puts "Number of files to update: #{file_count}"
  end

  private

  # Fetch the last_scan timestamp from the pipeline state using Polyphemus client
  def fetch_last_scan
    state = polyphemus_client.get_workflow_state(argo_id)
    state[:last_scan] ? Time.parse(state[:last_scan]) : Time.at(0) # Default to epoch if not found
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