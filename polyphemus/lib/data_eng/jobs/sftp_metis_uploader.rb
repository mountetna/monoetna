require_relative 'etl_job'
require_relative 'common'
require 'securerandom'

class SftpMetisUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithSlackNotifications
  include WithLogger

  METIS_FAILED_FILES_CSV = "metis_failed_files.csv"

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

  def name_regex
    Regexp.new("#{magic_string}(-|_)")
  end

  def bucket_name
    config['config']['bucket_name']
  end

  def metis_uid
    @metis_uid ||= SecureRandom.hex
  end

  def metis_root_path
    config['config']['metis_root_path']
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

  public

  def pre(context)
    run = polyphemus_client.get_run(project_name, run_id)

    raise "Run #{run_id} not found" unless run

    context[:files_to_update] = (run.dig("state","files_to_update") || []).map(&:symbolize_keys) 

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

    if context[:successful_files].any?
      msg += "Uploaded #{context[:successful_files].size} files:\n" +
        context[:successful_files].join("\n")
    end

    if context[:failed_files].any?
      msg += "Failed #{context[:failed_files].size} files:\n" +
        context[:failed_files].join("\n")

      polyphemus_client.update_run(project_name, run_id, {
        state: {
          metis_num_failed_files: context[:failed_files].size
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

