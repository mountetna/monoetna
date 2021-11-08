require_relative "../db_triage_file_notification_etl"
require_relative "../helpers"

class Polyphemus::SlackNotificationTriageBaseEtl < Polyphemus::DbTriageFileNotificationEtl
  include WithSlackNotifications

  private

  def serialize(records)
    records.map { |r| r[:name] }.join("\n")
  end
end
