require_relative 'etl_job'
require_relative 'sftp_config'
require_relative '../clients/sftp_client'

class SftpFileDiscoveryJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithLogger
  include WithSftpConfig

  private

  def restart_scan?
    !!DateTime.parse(initial_start_scan_time) rescue false
  end

  def sftp_client
    @sftp_client ||= SFTPClient.new(
      ingest_host,
      ingest_user,
      ingest_password
    )
  end

  public

  def pre(context)
    context[:start_time] = fetch_last_scan
    if override_interval
      context[:end_time] = context[:start_time] + override_interval
    else
      context[:end_time] = Time.now.to_i
    end
    true
  end

  def process(context)
    logger.info("Searching for files from #{context[:start_time]} to #{context[:end_time]}")
    files_to_update = sftp_client.search_files(
      ingest_root_path,
      file_regex,
      context[:start_time],
      context[:end_time ]
    )
    logger.info("Found #{files_to_update.size} files to update...")
    context[:files_to_update] = files_to_update
  end

  def post(context)
    polyphemus_client.update_run(project_name, run_id, {
      state: {
        num_files_to_update: context[:num_files_to_update],
        files_to_update: context[:files_to_update],
        start_time: context[:start_time],
        end_time: context[:end_time]
      }
    })
  end

  private

  def fetch_last_scan
    if restart_scan?
      logger.info("Restarting scan... at #{initial_start_scan_time}")

      return DateTime.parse(initial_start_scan_time).to_i
    end

    begin
      response = polyphemus_client.get_previous_state(
        project_name,
        workflow_config_id,
        state: [:end_time]
      )
      response["end_time"].to_i
    rescue Etna::Error => e
      logger.warn("Error fetching previous state")
      raise e
    end
  end
end
