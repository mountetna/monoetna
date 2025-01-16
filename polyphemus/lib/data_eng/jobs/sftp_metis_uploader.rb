require_relative 'etl_job'
require_relative 'common'
require_relative '../clients/sftp_client'

class SFTPMetisUploaderJob < Polyphemus::ETLJob
  include WithEtnaClients
  include WithLogger

  METIS_FAILED_FILES_CSV = "metis_failed_files.csv"

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
    @bucket_name = config['config']['bucket_name']
    @metis_uid = config['config']['metis_uid']
    @metis_root_path = config['config']['metis_root_path']
    @path_to_write_files = config['config']['path_to_write_files']
  end

  def pre(context)
    run = polyphemus_client.get_run(@project_name, run_id)
    unless run
        raise "Run #{run_id} not found"
    end
    context[:files_to_update] = CSV.foreach(run["state"]["files_to_update_path"], headers: true).map do |row|
      { path: row['path'], modified_time: row['modified_time'].to_i }
    end

    if context[:files_to_update].empty?
      logger.info("No new files to upload...")
      return false
    else
      return true
    end
  end

  def process(context)
    context[:failed_files] = []

    context[:files_to_update].each do |file|
      sftp_path = file[:path]
      modified_time = file[:modified_time]

      begin
        file_stream = @sftp_client.download_as_stream(sftp_path)
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
          project_name: @project_name,
          bucket_name: @bucket_name,
          metis_uid: @metis_uid,
        )
        metis_path = File.join(@metis_root_path, remove_dscolab_prefix(sftp_path))
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
      writable_dir = build_pipeline_state_dir(@path_to_write_files, run_id)
      context[:failed_files_path] = File.join(writable_dir, METIS_FAILED_FILES_CSV)

      logger.warn("Writing failed files to #{context[:failed_files_path]}")
      logger.warn("Found #{context[:failed_files].size} failed files")
      CSV.open(context[:failed_files_path], "wb") do |csv|
        csv << ["path"]
        context[:failed_files].each do |path|
        csv << [path]
        end
      end
      polyphemus_client.update_run(@project_name, @run_id, {
        state: {
          metis_num_failed_files: context[:failed_files].size,
          metis_failed_files_path: context[:failed_files_path],
        },
      })
    else
      logger.info("No failed files found")
    end
  end
end

