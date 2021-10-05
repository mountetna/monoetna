require_relative "../db_triage_file_notification_etl"
require_relative "./slack_notification_module"

class Polyphemus::SlackNotificationTriageBaseEtl < Polyphemus::DbTriageFileNotificationEtl
  include Polyphemus::WithSlackNotification

  private

  def serialize(records)
    records.map { |r| r[:name] }.join("\n")
  end
end
