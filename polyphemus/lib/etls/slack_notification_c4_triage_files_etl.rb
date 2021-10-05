require_relative "./slack_notification_triage_base_etl"

class Polyphemus::SlackNotificationC4TriageFilesEtl < Polyphemus::SlackNotificationTriageBaseEtl
  def initialize
    @project_name = "c4"
    @bucket_name = "triage"

    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      column_name: :triage_ingested_at,
    )
  end

  def process(cursor, records)
    slack(
      "C4 Triage from CAT ETL",
      "data-ingest-ping",
      "Downloaded #{records.length} files from the CAT to C4.\n#{serialize(records)}.\nThey are available at #{c4_config[:root]} and ready for triage."
    )
  end

  private

  def c4_config
    @c4_config ||= ssh_configs.find { |c| c[:alias] == "c4" }
  end

  def ssh_configs
    ingest_configs[:ssh]
  end

  def ingest_configs
    Polyphemus.instance.config(:ingest)
  end
end
