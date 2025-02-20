require_relative 'etl_job'
require_relative 'common'
require_relative '../clients/sftp_client'

class SftpDepositUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
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
    config['config']['deposit_root_path']
  end

  def sftp_client
    @sftp_client ||= SFTPClient.new(
      config["secrets"]["sftp_ingest_host"],
      config["secrets"]["sftp_ingest_user"],
      config["secrets"]["sftp_ingest_password"],
    )
  end

  def deposit_host
    config["secrets"]["sftp_deposit_host"]
  end

  def remote_ssh
    @remote_ssh ||= Etna::RemoteSSH.new(
      host: deposit_host,
      username: config["secrets"]["sftp_deposit_user"],
      password: config["secrets"]["sftp_deposit_password"],
      root: nil
    )
  end

  def pre(context)
    run = polyphemus_client.get_run(project_name, run_id)

    unless run
      raise "Run #{run_id} not found"
    end

    context[:files_to_update] = run["state"]["files_to_update"].map(&:symbolize_keys)

    if context[:files_to_update].empty?
      logger.info("No new files to upload...")
      return false
    end

    return true
  end

  def process(context)
    context[:failed_files] = []

    context[:files_to_update].each do |file|
      sftp_path = file[:path]
      modified_time = file[:modified_time]

      begin
        file_stream = sftp_client.download_as_stream(sftp_path)
        if file_stream.nil?
          logger.info("Failed to download #{sftp_path}")
            context[:failed_files] << sftp_path
          next
        end
      rescue StandardError => e
        logger.info("Failed to download from sftp server, #{sftp_path}: #{e.message}")
        next
      end

      begin
        deposit_path = File.join(deposit_root_path, sftp_path.gsub(name_regex,''))
        remote_ssh.file_upload(deposit_path, file_stream)
      rescue Etna::RemoteSSH::RemoteSSHError => e
        logger.warn("Failed to upload to deposit host #{deposit_host}: #{deposit_path}. Error: #{e.message}")
        context[:failed_files] << sftp_path 
      end
    end
  end

  def post(context)
    if context[:failed_files].any?
      logger.warn("Found #{context[:failed_files].size} failed files")

      context[:failed_files].each do |file_path|
        logger.warn(file_path)
      end

      polyphemus_client.update_run(project_name, run_id, {
        state: {
          deposit_num_failed_files: context[:failed_files].size
        },
      })
    else
      logger.info("No failed files found")
    end
  end
end
