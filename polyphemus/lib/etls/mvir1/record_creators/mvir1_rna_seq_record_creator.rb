require_relative "../../../helpers"

class Polyphemus::MvirRnaSeqRecordCreator
  include WithLogger
  include WithEtnaClients
  include WithSlackNotifications

  PATH_REGEX = /^bulk_RNASeq\/raw\/(?<record_name>.*)/
  TIMEPOINT_REGEX = /^(?<timepoint>MVIR1-HS\d+-DN?\d+)[A-Z]+.*/

  def initialize
  end

  def process(cursor, folders)
    create_records(cursor, folders)

    logger.info("Done")
  end

  def create_records(cursor, folders)
    folders.each do |folder|
      record_name = ::File.basename(folder.folder_path)

      attrs = {
        timepoint: timepoint(record_name),
      }

      if attrs[:timepoint].nil?
        notify_slack(
          "Skipping bulk rna_seq record without valid timepoint #{record_name}.",
          channel: "data-ingest-errors",
        )
        logger.info("No timepoint, skipping #{record_name}.")
        next
      end

      begin
        update_request = Etna::Clients::Magma::UpdateRequest.new(
          project_name: cursor[:project_name],
        )
        update_request.update_revision("rna_seq", record_name, attrs)
        logger.info("Creating rna_seq records: #{update_request.revisions["rna_seq"].keys.join(",")}")
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
end
