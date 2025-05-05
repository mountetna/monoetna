require_relative 'etl_job'
require_relative 'sftp_config'
require 'securerandom'

class SftpMetisUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithSlackNotifications
  include WithLogger
  include WithSftpConfig

  METIS_FAILED_FILES_CSV = "metis_failed_files.csv"

  private

  def metis_uid
    @metis_uid ||= SecureRandom.hex
  end

  public

  def pre(context)
    run = polyphemus_client.get_run(project_name, run_id)

    raise "Run #{run_id} not found" unless run

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

  def metis_filesystem
    Etna::Filesystem::Metis.new(
      project_name: project_name,
      bucket_name: bucket_name,
      metis_client: metis_client,
      root: metis_root_path
    )
  end

  def ingest_filesystem
    Etna::Filesystem::SftpFilesystem.new(
      username: config["secrets"]["sftp_ingest_user"],
      password: config["secrets"]["sftp_ingest_password"],
      host: ingest_host
    )
  end

  def process(context)
    context[:failed_files] = []
    context[:successful_files] = []

    workflow = Etna::Clients::Metis::IngestMetisDataWorkflow.new(
      metis_filesystem: metis_filesystem,
      ingest_filesystem: ingest_filesystem,
      logger: nil,
    )

    logger.info("Uploading #{context[:files_to_update].size} to metis at #{metis_root_path}")

    workflow.copy_files(context[:files_to_update].map do |file|
      [ file[:path], file[:path].gsub(name_regex, '') ]
    end) do |filename, success|
      if success
        context[:successful_files] << filename
      else
        logger.warn("Failed to upload to metis: #{filename}")
        context[:failed_files] << filename
      end
    end
  end

  def post(context)
    msg =
      "Finished uploading #{context[:files_to_update].size} files " +
      "from #{ingest_host} to Metis for #{project_name}. " +
      "Please check #{metis_root_path} path in bucket #{bucket_name}."

    msg += "\n" 

    state = {}

    if context[:successful_files]&.any?
      msg += "Uploaded #{context[:successful_files].size} files:\n" +
        context[:successful_files].join("\n")

      state[:metis_successful_files] = context[:successful_files]
    end

    if context[:failed_files]&.any?
      msg += "Failed #{context[:failed_files].size} files:\n" +
        context[:failed_files].join("\n")

      logger.warn("Found #{context[:failed_files].size} failed files")

      context[:failed_files].each do |file_path|
        logger.warn(file_path)
      end

      state[:metis_num_failed_files] = context[:failed_files].size
    end
    
    unless state.empty?
      polyphemus_client.update_run(project_name, run_id, { state: state })
      notify_slack(msg, channel: notification_channel, webhook_url: notification_webhook_url)
    end
  end

  private

  def fetch_successful_files
    begin
      response = polyphemus_client.get_previous_state(
        project_name,
        workflow_config_id,
        state: [:metis_successful_files],
        collect: true
      )

      files = (response["metis_successful_files"] || []).flatten.map do |filename|
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

