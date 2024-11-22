require_relative 'etl_job'
require_relative 'common'
require_relative '../clients/sftp_client'

class SFTPFileDiscoveryJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithLogger

  SFTP_FILES_TO_UPDATE_CSV = "sftp_files_to_update.csv"

  def initialize(config, runtime_config)
    @project_name = config['project_name']
    @sftp_client = SFTPClient.new(
      config["secrets"]["sftp_host"],
      config["secrets"]["sftp_user"],
      config["secrets"]["sftp_password"],
    )

    # Mandatory params
    @workflow_config_id = config['config_id']
    @workflow_version = config['version_number']
    @file_regex = config['config']['file_regex']
    @sftp_root_dir = config['config']['sftp_root_dir']
    @initial_start_scan_time = config['config']['initial_start_scan_time']
    @path_to_write_files = config['config']['path_to_write_files']

    # Optional params
    @interval = config['config']['interval'] || nil
  end

  def pre(context)
    context[:start_time] = fetch_last_scan
    if @interval
      context[:end_time] = context[:start_time] + @interval
    else
      context[:end_time] = Time.now.to_i
    end
    true
  end

  def process(context)
    logger.info("Searching for files from #{context[:start_time]} to #{context[:end_time]}")
    files_to_update = @sftp_client.search_files(@sftp_root_dir, @file_regex, context[:start_time], context[:end_time ])
    if files_to_update.empty?
      context[:num_files_to_update] = 0
      logger.info("No files found to update...")
    else
      context[:num_files_to_update] = files_to_update.size
      writable_dir = build_pipeline_state_dir(@path_to_write_files, run_id)
      context[:files_to_update_path] = write_csv(writable_dir, files_to_update)
      logger.info("Found #{files_to_update.size} files to update...")
    end
  end

  def post(context)
    polyphemus_client.update_run(@project_name, @run_id, @workflow_config_id, @workflow_version, {
      state: {
        num_files_to_update: context[:num_files_to_update],
        files_to_update_path: context[:files_to_update_path],
        start_time: context[:start_time],
        end_time: context[:end_time]
      },
    })
  end

  private

  # Fetch the last_scan timestamp from the pipeline state using Polyphemus client
  def fetch_last_scan
    run = polyphemus_client.get_previous_run(@project_name, @workflow_config_id, @workflow_version)
    # No previous run exists, this is the first run ever
    if run.empty?
      logger.info("No previous run exists, this is the first run ever...")
      @initial_start_scan_time
    else
      logger.info("Previous run exists, using end_time from previous run...")
      run["state"]["end_time"]
    end
  end

  def write_csv(dir_path, files_to_update)
    filepath = File.join(dir_path, SFTP_FILES_TO_UPDATE_CSV)
    CSV.open(filepath, "wb") do |csv|
      csv << ["path", "modified_time"]
      files_to_update.each do |file|
        csv << [file[:path], file[:modified_time]]
      end
    end
    filepath
  end
end
