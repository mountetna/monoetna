require_relative "../db_triage_file_etl"

class Polyphemus::SftpIngestCatToC4TriageFilesEtl < Polyphemus::DbTriageFileEtl
  def initialize
    @project_name = "c4"
    @bucket_name = "triage"

    raise "No CAT configuration provided" unless cat_config
    raise "No C4 configuration provided" unless c4_config

    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      limit: 20,
      column_name: :triage_ingested_at,
    )
  end

  def process(cursor, records)
    cat_records = records.select do |record|
      record[:host] == cat_config[:host]
    end

    logger.info("Ingesting files from #{cat_config[:host]}: #{cat_records.map { |c| c[:name] }.join(", ")}...")

    cat_records.each do |record|
      c4_connection.lftp_get(
        host: cat_config[:host],
        username: cat_config[:username],
        password: cat_config[:password],
        remote_filename: record[:name],
      ) do |file_name|
        logger.info("#{file_name} completed.")
        update_ingested_timestamp(file_name)
        `/bin/post-to-slack "C4 Triage from CAT ETL" "data-ingest-ping" "Successfully downloaded #{file_name} from the CAT to C4. It is available at #{c4_config[:root]} and ready for triage." || true`
      end
    end

    logger.info("Done")
  end

  private

  def cat_config
    @cat_config ||= sftp_configs.find { |c| c[:alias] == "cat" }
  end

  def c4_config
    @c4_config ||= ssh_configs.find { |c| c[:alias] == "c4" }
  end

  def sftp_configs
    ingest_configs[:sftp]
  end

  def ssh_configs
    ingest_configs[:ssh]
  end

  def ingest_configs
    Polyphemus.instance.config(:ingest)
  end

  def update_ingested_timestamp(file_name)
    Polyphemus::IngestFile.where(
      host: cat_config[:host],
      name: file_name,
    ).first.update(triage_ingested_at: DateTime.now)
  end

  def c4_connection
    @c4 ||= Etna::RemoteSSH.new(
      host: c4_config[:host],
      username: c4_config[:username],
      password: c4_config[:password],
      root: c4_config[:root],
    )
  end
end
