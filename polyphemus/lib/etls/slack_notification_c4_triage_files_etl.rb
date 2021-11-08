require_relative "./slack_notification_triage_base_etl"
require_relative "./ingest_config_module"

class Polyphemus::SlackNotificationC4TriageFilesEtl < Polyphemus::SlackNotificationTriageBaseEtl
  include Polyphemus::WithIngestConfig

  def initialize
    @project_name = "c4"
    @bucket_name = "triage"

    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      column_name: :triage_ingested_at,
    )
  end

  def process(cursor, records)
    notify_slack(
      "Downloaded #{records.length} files from the CAT to C4.\n#{serialize(records)}.\nThey are available at #{c4_config[:root]} and ready for triage.",
      channel: "data-ingest-ping",
    )
  end
end
