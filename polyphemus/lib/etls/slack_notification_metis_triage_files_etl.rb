require_relative "./slack_notification_triage_base_etl"

class Polyphemus::SlackNotificationMetisTriageFilesEtl < Polyphemus::SlackNotificationTriageBaseEtl
  def initialize
    @project_name = "triage"
    @bucket_name = "waiting_room"

    super(
      project_bucket_pairs: [[@project_name, @bucket_name]],
      column_name: :archive_ingested_at,
    )
  end

  def process(cursor, records)
    slack(
      "Metis Triage from CAT ETL",
      "data-ingest-ping",
      "Downloaded #{records.length} files from the CAT to Metis.\n#{serialize(records)}.\nThey are available in #{@project_name} / #{@bucket_name}."
    )
  end
end
