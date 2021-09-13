require_relative "../db_triage_file_etl"

class Polyphemus::SftpIngestMetisTriageFilesEtl < Polyphemus::DbTriageFileEtl
  def initialize
    @project_name = "triage"
    @bucket_name = "waiting_room"
    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      limit: 20,
      column_name: :archive_ingested_at,
    )
  end

  def process(cursor, records)
    records_by_host = records.group_by do |record|
      record[:host]
    end

    sftp_configs.each do |conf|
      host = conf[:host]
      files_for_host = records_by_host[host]
      file_names = files_for_host.map { |file| file[:name] }

      logger.info("Ingesting #{file_names.length} files from #{host}: #{file_names.join(", ")}...")

      workflow = Etna::Clients::Metis::IngestMetisDataWorkflow.new(
        metis_filesystem: metis_filesystem(conf),
        ingest_filesystem: ingest_filesystem(conf),
        logger: logger,
      )
      workflow.copy_files(file_names) do |file_name|
        puts "#{file_name} finished uploading."
        update_ingested_timestamp(host, file_name)
        `/bin/post-to-slack.sh "Metis Archive from CAT ETL" "data-ingest-ping" "Successfully archived #{file_name} from the CAT on Metis." || true`
      end
    end

    logger.info("Done")
  end

  private

  def sftp_configs
    Polyphemus.instance.config(:ingest)[:sftp]
  end

  def update_ingested_timestamp(host, file_name)
    Polyphemus::IngestFile.where(host: host, name: file_name).first.update(
      archive_ingested_at: DateTime.now,
    )
  end

  def ingest_filesystem(configuration)
    Etna::Filesystem::SftpFilesystem.new(**configuration)
  end

  def metis_filesystem(configuration)
    # Set up the root to be the configured host, so that
    #   files will be in a name-spaced directory on Metis.
    Etna::Filesystem::Metis.new(
      metis_client: metis_client,
      project_name: @project_name,
      bucket_name: @bucket_name,
      root: configuration[:host],
    )
  end
end
