require_relative "../../../helpers"
require_relative "../../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqAndPlateRecordCreator
  include WithLogger
  include WithEtnaClients
  include WithSlackNotifications

  PATH_REGEX = /.*\/(?<plate>plate\d+)_.*\/output\/(?<record_name>.*)/
  SAMPLE_NAME_REGEX = /^(?<sample_name>IPI.*\.[A-Z]+\d)\..*/

  def initialize
    @helper = IpiHelper.new
  end

  def process(cursor, folders)
    # We can extract the plate# from the folder_path:
    #   bulkRNASeq/processed/plate1_rnaseq_new/output/blahblahrecordname
    # We will also have to format the Control record names to match IPI validation.

    @plate_names = plate_names(folders)

    logger.info("Ensuring plates exist: #{@plate_names.join(",")}")
    ensure_plates(cursor, @plate_names)
    create_records(cursor, folders)

    logger.info("Done")
  end

  def plate_names(folders)
    folders.map do |folder|
      plate(folder.folder_path)
    end.uniq
  end

  def plate(folder_path)
    match = folder_path.match(PATH_REGEX)

    match[:plate].capitalize
  end

  def ensure_plates(cursor, plate_names)
    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: cursor[:project_name],
    )
    plate_names.each do |plate_name|
      update_request.update_revision("rna_seq_plate", plate_name, {
        "project" => "UCSF Immunoprofiler",
      })
    end
    magma_client.update_json(update_request)
  end

  def create_records(cursor, folders)
    folders.each do |folder|
      plate_name = plate(folder.folder_path)
      folder_name = ::File.basename(folder.folder_path)

      next if @helper.is_non_cancer_sample?(folder_name)

      record_name = @helper.is_control?(folder_name) ?
        @helper.control_name(folder_name) :
        @helper.corrected_rna_seq_tube_name(folder_name)

      attrs = {
        rna_seq_plate: plate_name,
      }

      unless @helper.is_control?(folder_name)
        attrs[:sample] = sample_name(record_name)
        if attrs[:sample].nil?
          notify_slack(
            "Skipping non control record without valid sample name #{record_name}.",
            channel: 'data-ingest-errors')
          next
        end
      end

      begin
        update_request = Etna::Clients::Magma::UpdateRequest.new(
          project_name: cursor[:project_name],
        )
        update_request.update_revision("rna_seq", record_name, attrs)
        logger.info("Creating rna_seq records: #{update_request.revisions["rna_seq"].keys.join(",")}")
        magma_client.update_json(update_request)
      rescue Exception => e
        notify_slack("Error creating IPI rna_seq record #{record_name}.\n#{e.message}.", channel: 'data-ingest-errors')
        logger.log_error(e)
      end
    end
  end

  def sample_name(rna_seq_record_name)
    match = rna_seq_record_name.match(SAMPLE_NAME_REGEX)
    match ? match[:sample_name] : nil
  end
end
