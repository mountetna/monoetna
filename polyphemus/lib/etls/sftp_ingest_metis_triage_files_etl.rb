require_relative "../db_triage_file_etl"

class Polyphemus::SftpIngestMetisTriageFilesEtl < Polyphemus::DbTriageFileEtl
  def initialize
    @project_name = "triage"
    @bucket_name = "waiting_room"
    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      limit: 5,
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

      logger.info("Ingesting files from #{host}: #{file_names.join(", ")}...")

      workflow = Etna::Clients::Metis::IngestMetisDataWorkflow.new(
        metis_filesystem: metis_filesystem(conf),
        ingest_filesystem: ingest_filesystem(conf),
        logger: logger,
      )
      workflow.copy_files(file_names)

      update_ingested_timestamp(files_for_host)
    end

    logger.info("Done")
  end

  private

  def sftp_configs
    Polyphemus.instance.config(:ingest)[:sftp]
  end

  def update_ingested_timestamp(file_records)
    Polyphemus::IngestFile.where(id: file_records.map { |f| f[:id] })
      .all do |file|
      file.update(ingested_at: DateTime.now)
    end
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
