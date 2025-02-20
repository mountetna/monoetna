require_relative 'etl_job'
require_relative 'common'
require_relative '../clients/sftp_client'

class SftpFileDiscoveryJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithLogger

  private

  def project_name
    config['project_name']
  end

  def workflow_config_id
    config['config_id']
  end

  def workflow_version
    config['version_number']
  end

  def magic_string
    config['config']['magic_string']
  end

  def file_regex
    Regexp.new("#{magic_string}(-|_).*")
  end

  def ingest_root_path
    config['config']['ingest_root_path']
  end

  def restart_scan?
    !!runtime_config['config']['restart_scan']
  end

  def initial_start_scan_time
    runtime_config['config']['initial_start_scan_time']
  end

  def override_interval
    runtime_config['config']['override_interval']
  end

  def sftp_client
    @sftp_client ||= SFTPClient.new(
      config["secrets"]["sftp_ingest_host"],
      config["secrets"]["sftp_ingest_user"],
      config["secrets"]["sftp_ingest_password"],
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
      logger.info("Restarting scan... at #{Time.at(initial_start_scan_time).strftime('%Y-%m-%d %H:%M:%S')}")

      return initial_start_scan_time
    else
      response = polyphemus_client.get_previous_state(
       project_name,
       workflow_config_id,
       workflow_version,
       state: [:end_time]
      )
      if response[:error]
        logger.warn("Error fetching previous state: #{response[:error]}")
        raise StandardError, response[:error]
      end
      response["end_time"]
    end
  end
end
