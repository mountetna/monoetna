# Insert and remove file pointers to / from the ingestion table
#    based on a scan of an SFTP server.
require_relative "../rsync_files_etl"

class Polyphemus::SyncCatSftpFilesEtl < Polyphemus::RsyncFilesEtl
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

  def host
    sftp_conf[:host]
  end

  def sftp_conf
    conf = sftp_configs.find do |conf|
      "cat" == conf[:alias].downcase
    end

    raise "No configuration found for the CAT. Make sure an `:ingest -> :sftp -> :alias: cat` entry exists in config.yml." if conf.nil?

    conf
  end

  def sftp_configs
    Polyphemus.instance.config(:ingest)[:sftp]
  end

  def update_ingest_files(records)
    # First, we filter out any files that have already been ingested
    records = uningested_records(records)

    # Next, we add new files that do not have a record in Polyphemus::IngestFiles
    add_new_ingest_files(records)

    # Finally, we'll remove any IngestFile records that no longer exist
    #   on the CAT.
    remove_deleted_ingest_files(records)
  end

  def ingest_files_for_host_query
    Polyphemus::IngestFile.where(host: host)
  end

  def uningested_records(records)
    records.select do |record|
      !ingested_file_names.include?(record.filename)
    end
  end

  def ingested_file_names
    @ingested_file_names ||= ingest_files_for_host_query.exclude(ingested_at: nil).all.map do |file|
      file[:name]
    end
  end

  def add_new_ingest_files(records)
    existing_names = existing_file_names(records)

    new_records = records.select do |record|
      !existing_names.include?(record.filename)
    end

    new_records.each do |record|
      Polyphemus::IngestFile.create(
        name: record.filename,
        host: host,
        should_ingest: false,
      )
    end
  end

  def existing_file_names(records)
    ingest_files_for_host_query.where(
      name: records.map { |r| r.filename },
    ).all.map { |f| f[:name] }
  end

  def remove_deleted_ingest_files(records)
    current_cat_filenames = records.map do |record|
      record.filename
    end

    ingest_files_for_host_query.where(ingested_at: nil)
      .exclude(name: current_cat_filenames)
      .all do |file|
      file.delete
    end
  end
end
