require_relative 'etl_job'
require_relative 'common'
require_relative '../clients/sftp_client'
require 'securerandom'

class SftpMetisUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
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

  def sftp_client
    @sftp_client ||= SFTPClient.new(
      config["secrets"]["sftp_ingest_host"],
      config["secrets"]["sftp_ingest_user"],
      config["secrets"]["sftp_ingest_password"],
    )
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

  def process(context)
    context[:failed_files] = []

    context[:files_to_update].each do |file|
      sftp_path = file[:path]
      modified_time = file[:modified_time]

      begin
        file_stream = sftp_client.download_as_stream(sftp_path)
        if file_stream.nil?
          logger.info("Failed to download #{sftp_path}")
          next
        end
      rescue StandardError => e
        logger.info("Failed to download from sftp server, #{sftp_path}: #{e.message}")
        context[:failed_files] << sftp_path
        next
      end

      begin
        uploader = Etna::Clients::Metis::MetisUploadWorkflow.new(
          metis_client: metis_client,
          project_name: project_name,
          bucket_name: bucket_name,
          metis_uid: metis_uid,
        )

        metis_path = File.join(metis_root_path, sftp_path.gsub(name_regex,''))
        uploader.do_upload(
          Etna::Clients::Metis::MetisUploadWorkflow::StreamingIOUpload.new(
            readable_io: file_stream,
            size_hint: file_stream.size,
          ),
          metis_path
        )
      rescue StandardError => e
        logger.warn("Failed to upload to metis: #{metis_path}. Error: #{e.message}")
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
          metis_num_failed_files: context[:failed_files].size
        },
      })
    else
      logger.info("No failed files found")
    end
  end
end

