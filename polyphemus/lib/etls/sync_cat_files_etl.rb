# Insert and remove file pointers to / from the ingestion table
#    based on an rsync scan of the CAT SFTP server.
require_relative "../rsync_files_etl"

class Polyphemus::SyncCatFilesEtl < Polyphemus::RsyncFilesEtl
  def initialize
    super(
      project_bucket_pairs: [["triage", "waiting_room"]],
      host: sftp_conf[:host],
      username: sftp_conf[:username],
      root: sftp_conf[:root],
      password: sftp_conf[:password],
    )
  end

  def process(cursor, records)
    logger.info("Processing #{records.length} files from CAT.")

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

  def update_ingest_files(all_records)
    # First, we filter out any files that have already been ingested.
    # This way they don't get re-ingested, and we'll also have
    #   a log of when stuff made it to Metis.
    uningested_records = uningested_records(all_records)

    logger.info("Found #{uningested_records.length} uningested files.")

    # Next, we add new files matching our magic string
    #  that do not have a record in Polyphemus::IngestFiles
    add_new_ingest_files(uningested_records)

    # Finally, we'll remove any IngestFile records that no longer exist
    #   on the CAT.
    remove_deleted_ingest_files(all_records)
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
    @ingested_file_names ||= ingest_files_for_host_query.exclude(triage_ingested_at: nil, archive_ingested_at: nil)
      .all.map do |file|
      file[:name]
    end
  end

  def add_new_ingest_files(records)
    existing_names = existing_file_names(records)

    new_records = records.select do |record|
      !existing_names.include?(record.filename) && record.filename =~ Regexp.new(filepath_regex)
    end

    logger.info("Found #{new_records.length} new files matching the regex: #{filepath_regex}.")

    Polyphemus::IngestFile.multi_insert(new_records.map do |record|
      {
        created_at: DateTime.now,
        updated_at: DateTime.now,
        name: record.filename,
        host: host,
        should_ingest: true,
      }
    end)
  end

  def filepath_regex
    Polyphemus.instance.config(:ingest)[:filepath_regex] || ""
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

    files_to_remove_query = ingest_files_for_host_query.exclude(name: current_cat_filenames)

    logger.info("Found #{files_to_remove_query.count} matching files removed from the CAT.")

    files_to_remove_query.all do |file|
      file.update(removed_from_source: true)
    end
  end
end
