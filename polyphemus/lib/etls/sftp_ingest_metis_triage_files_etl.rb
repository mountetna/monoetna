require_relative "../db_triage_file_etl"

class Polyphemus::SftpIngestMetisTriageFilesEtl < Polyphemus::DbTriageFileEtl
  def initialize
    @project_name = "triage"
    @bucket_name = "waiting_room"
    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      limit: 20,
      table_name: "metis_ingest_triage",
    )
  end

  def process(cursor, records)
    workflow = Etna::Clients::Metis::IngestMetisDataWorkflow.new(
      metis_filesystem: metis_filesystem,
      ingest_filesystem: ingest_filesystem,
      logger: logger,
    )

    host_matches = records.select do |record|
      host == record[:host]
    end

    logger.warn("#{records.length - host_matches.length} triage file(s) are from a different SFTP host than the current configuration: #{host}") if host_matches.length != records.length

    logger.info("Ingesting files from #{host}: #{host_matches.map { |r| r[:name] }.join(", ")}...")

    workflow.ingest_files(host_matches.map { |m| m[:name] })

    logger.info("Done")
  end

  private

  def config
    Polyphemus.instance.config(:ingest)[:sftp]
  end

  def host
    config[:host]
  end

  def ingest_filesystem
    Etna::Filesystem::SftpFilesystem.new(
      host: host,
      username: config[:username],
      password: config[:password],
    )
  end

  def metis_filesystem
    Etna::Filesystem::Metis.new(
      metis_client: metis_client,
      project_name: @project_name,
      bucket_name: @bucket_name,
    )
  end
end
