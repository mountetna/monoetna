# Insert and remove file pointers to / from the ingestion table
#    based on a scan of an SFTP server.
require_relative "../rsync_ingest_file_etl"

class Polyphemus::SyncCatSftpFilesEtl < Polyphemus::RsyncIngestFileEtl
  def initialize
    @project_name = "triage"
    @bucket_name = "waiting_room"
    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      host: sftp_conf[:host],
      username: sftp_conf[:username],
      root: sftp_conf[:root],
      password: sftp_conf[:password],
    )
  end

  def process(cursor, records)
    logger.info("Ingesting files from CAT: #{records.map { |r| r.filename }.join(", ")}...")

    update_ingest_files(records)

    logger.info("Done")
  end

  private

  def sftp_conf
    conf = sftp_configs.find do |conf|
      "cat" == conf[:alias]
    end

    raise "No configuration found for the CAT." if conf.nil?

    conf
  end

  def sftp_configs
    Polyphemus.instance.config(:ingest)[:sftp]
  end

  def update_ingest_files(records)
    binding.pry
    # First, we filter out any files that have already been ingested
    records = uningested_records(records)

    # Next, we add new files that do not have a record in IngestFiles

    # Finally, we'll remove any IngestFile records that no longer exist
    #   on the CAT.

  end

  def uningested_records(records)
    ingested_file_names = Polyphemus::IngestFile.exclude(ingested_at: nil).all.map do |file|
    end
  end
end
