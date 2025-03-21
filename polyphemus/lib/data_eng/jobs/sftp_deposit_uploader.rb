require_relative 'etl_job'
require_relative 'sftp_config'

class SftpDepositUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithSlackNotifications
  include WithLogger
  include WithSftpConfig

  DEPOSIT_FAILED_FILES_CSV = "deposit_failed_files.csv"

  def remote_ssh
    @remote_ssh ||= Etna::RemoteSSH.new(
      host: deposit_host,
      username: deposit_user,
      password: deposit_password,
      root: deposit_root_path
    )
  end

  def pre(context)
    run = polyphemus_client.get_run(project_name, run_id)

    unless run
      raise "Run #{run_id} not found"
    end

    context[:files_to_ignore] = fetch_successful_files

    context[:files_to_update] = (run.dig("state","files_to_update") || []).map(&:symbolize_keys).reject do |file|
      context[:files_to_ignore].include?(file[:path])
    end

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
          username: ingest_user,
          password: ingest_password,
          host: ingest_host,
          remote_filename: sftp_path,
          local_filename: deposit_path
        )
        context[:successful_files] << sftp_path
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

    state = {}

    if context[:successful_files]&.any?
      msg += "Uploaded #{context[:successful_files].size} files:\n" +
        context[:successful_files].join("\n")

      state[:deposit_successful_files] = context[:successful_files]
    end

    if context[:failed_files]&.any?
      msg += "Failed #{context[:failed_files].size} files:\n" +
        context[:failed_files].join("\n")

      logger.warn("Found #{context[:failed_files].size} failed files")

      context[:failed_files].each do |file_path|
        logger.warn(file_path)
      end

      state[:deposit_num_failed_files] = context[:failed_files].size
    end
    
    unless state.empty?
      polyphemus_client.update_run(project_name, run_id, { state: state })
    end

    notify_slack(msg, channel: notification_channel, webhook_url: notification_webhook_url)
  end

  private

  def fetch_successful_files
    begin
      response = polyphemus_client.get_previous_state(
        project_name,
        workflow_config_id,
        state: [:deposit_successful_files],
        collect: true
      )
      files = (response["deposit_successful_files"] || []).flatten.map do |filename|
        override_root_path ? [
          filename,
          filename.sub(/^#{raw_ingest_root_path}/, override_root_path)
        ] : filename
      end.flatten

      return Set.new(files)
    rescue Etna::Error => e
      logger.warn("Error fetching previous state")
      raise e
    end
  end
end
