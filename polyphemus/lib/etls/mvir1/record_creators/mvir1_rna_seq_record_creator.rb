require_relative "../../../helpers"

class Polyphemus::MvirRnaSeqRecordCreator
  include WithLogger
  include WithEtnaClients
  include WithSlackNotifications

  PATH_REGEX = /^bulk_RNASeq\/raw\/(?<record_name>.*)/
  TIMEPOINT_REGEX = /^(?<timepoint>MVIR1-HS\d+-DN?\d+)[A-Z]+.*/
  PATIENT_REGEX = /^(?<patient>MVIR1-HS\d+)-D.*/
  TIMEPOINT_DAY_REGEX = /.*-DN?(?<day>\d+)$/

  def initialize
  end

  def process(cursor, folders)
    create_records(cursor, folders)

    logger.info("Done")
  end

  def create_records(cursor, folders)
    folders.each do |folder|
      record_name = ::File.basename(folder.folder_path)

      timepoint_id = timepoint(record_name)
      patient_id = patient(timepoint_id)

      attrs = {
        timepoint: timepoint_id,
      }

      timepoint_attrs = {
        patient: patient_id,
        day: day(timepoint_id)
      }

      if timepoint_id.nil?
        notify_slack(
          "Skipping bulk rna_seq record without valid timepoint #{record_name}.",
          channel: "data-ingest-errors",
        )
        logger.info("No timepoint, skipping #{record_name}.")
        next
      end

      if patient_id.nil?
        notify_slack(
          "Skipping bulk rna_seq record without valid patient #{record_name}.",
          channel: "data-ingest-errors",
        )
        logger.info("No patient, skipping #{record_name}.")
        next
      end

      begin
        update_request = Etna::Clients::Magma::UpdateRequest.new(
          project_name: cursor[:project_name],
        )
        update_request.update_revision("rna_seq", record_name, attrs)
        update_request.update_revision("timepoint", timepoint_id, timepoint_attrs)
        logger.info("Creating rna_seq records: #{update_request.revisions["rna_seq"].keys.join(",")}")
        logger.info("Creating timepoint records: #{update_request.revisions["timepoint"].keys.join(",")}")
        magma_client.update_json(update_request)
      rescue Exception => e
        notify_slack("Error creating MVIR1 rna_seq record #{record_name}.\n#{e.message}.", channel: "data-ingest-errors")
        logger.log_error(e)
      end
    end
  end

  def timepoint(rna_seq_record_name)
    match = rna_seq_record_name.match(TIMEPOINT_REGEX)
    match ? match[:timepoint] : nil
  end

  def patient(timepoint)
    match = timepoint.match(PATIENT_REGEX)
    match ? match[:patient] : nil
  end

  def day(timepoint)
    match = timepoint.match(TIMEPOINT_DAY_REGEX)

    return nil unless match

    is_negative = timepoint.match(/DN/)

    return match[:day].to_i unless is_negative

    -(match[:day].to_i)
  end
end
