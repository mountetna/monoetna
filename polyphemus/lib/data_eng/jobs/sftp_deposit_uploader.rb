require_relative 'etl_job'
require_relative 'common'

class SftpDepositUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithSlackNotifications
  include WithLogger

  DEPOSIT_FAILED_FILES_CSV = "deposit_failed_files.csv"

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

  def name_regex
    Regexp.new("#{magic_string}(-|_)")
  end

  def deposit_root_path
    # Mandatory params
    ::File.join('/', config['config']['deposit_root_path'] || '')
  end

  def deposit_host
    config["secrets"]["sftp_deposit_host"]
  end

  def ingest_host
    config["secrets"]["sftp_ingest_host"]
  end

  def ingest_host
    config["secrets"]["sftp_ingest_host"]
  end

  def notification_channel
    config["config"]["notification_channel"]
  end

  def notification_webhook_url
    config["secrets"]["notification_webhook_url"]
  end

  def remote_ssh
    @remote_ssh ||= Etna::RemoteSSH.new(
      host: deposit_host,
      username: config["secrets"]["sftp_deposit_user"],
      password: config["secrets"]["sftp_deposit_password"],
      root: deposit_root_path
    )
  end

  def pre(context)
    run = polyphemus_client.get_run(project_name, run_id)

    unless run
      raise "Run #{run_id} not found"
    end

    context[:files_to_update] = (run.dig("state","files_to_update") || []).map(&:symbolize_keys) 

    if context[:files_to_update].empty?
      logger.info("No new files to upload...")
      return false
    end

    return true
  end

  def process(context)
    context[:failed_files] = []
    context[:successful_files] = []

    logger.info("Uploading #{context[:files_to_update].size} to #{deposit_host} at #{deposit_root_path}")

    context[:files_to_update].each do |file|
      sftp_path = file[:path]
      modified_time = file[:modified_time]
      deposit_path = ::File.join(deposit_root_path, ingest_host, sftp_path.gsub(name_regex,''))

      begin
        remote_ssh.lftp_get(
          username: config["secrets"]["sftp_ingest_user"],
          password: config["secrets"]["sftp_ingest_password"],
          host: ingest_host,
          remote_filename: sftp_path,
          local_filename: deposit_path
        )
        context[:successful_files] << deposit_path
      rescue Etna::RemoteSSH::RemoteSSHError => e
        logger.warn("Failed to upload to deposit host #{deposit_host}: #{sftp_path}. Error: #{e.message}")
        context[:failed_files] << sftp_path 
      end
    end
  end

  def post(context)
    msg =
      "Finished uploading #{context[:files_to_update].size} files " +
      "from #{ingest_host} to #{deposit_host} for #{project_name}. " +
      "Please check #{deposit_root_path} path."

    msg += "\n" 

    if context[:successful_files].any?
      msg += "Uploaded #{context[:successful_files].size} files:\n" +
        context[:successful_files].join("\n")
    end

    if context[:failed_files].any?
      msg += "Failed #{context[:failed_files].size} files:\n" +
        context[:failed_files].join("\n")

      polyphemus_client.update_run(project_name, run_id, {
        state: {
          deposit_num_failed_files: context[:failed_files].size
        },
      })

      logger.warn("Found #{context[:failed_files].size} failed files")

      context[:failed_files].each do |file_path|
        logger.warn(file_path)
      end
    end
    
    notify_slack(msg, channel: notification_channel, webhook_url: notification_webhook_url)
  end
end
