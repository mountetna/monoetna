require_relative "../../../helpers"
require_relative "../../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqAndPlateRecordCreator
  include WithLogger
  include WithEtnaClients

  PATH_REGEX = /.*\/(?<plate>plate\d+)_.*\/output\/(?<record_name>.*)/
  SAMPLE_NAME_REGEX = /^(?<sample_name>IPI.*\.[A-Z]+\d)\..*/
  PROJECT = "ipi"
  BUCKET = "data"

  def initialize
    @helper = IpiHelper.new
  end

  def create(folder_path)
    # We can extract the plate# from the folder_path:
    #   bulkRNASeq/processed/plate1_rnaseq_new/output/blahblahrecordname
    # We will also have to format the Control record names to match IPI validation.

    folder_name = ::File.basename(folder_path)
    return if @helper.is_non_cancer_sample?(folder_name)

    plate_name = plate(folder_path)

    logger.info("Ensuring plate exists: #{plate_name}")
    ensure_plate(plate_name)
    create_record(folder_path)

    logger.info("Done")
  end

  def plate(folder_path)
    match = folder_path.match(PATH_REGEX)

    match[:plate].capitalize
  end

  def ensure_plate(plate_name)
    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )
    update_request.update_revision("rna_seq_plate", plate_name, {
      "project" => "UCSF Immunoprofiler",
    })
    magma_client.update_json(update_request)
  end

  def create_record(folder_path)
    plate_name = plate(folder_path)
    folder_name = ::File.basename(folder_path)

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )

    record_name = @helper.is_control?(folder_name) ?
      @helper.control_name(folder_name) :
      @helper.corrected_rna_seq_tube_name(folder_name)

    attrs = {
      rna_seq_plate: plate_name,
    }

    attrs[:sample] = sample_name(record_name) unless @helper.is_control?(folder_name)

    update_request.update_revision("rna_seq", record_name, attrs)

    logger.info("Creating rna_seq records: #{update_request.revisions["rna_seq"].keys.join(",")}")
    magma_client.update_json(update_request)
  end

  def sample_name(rna_seq_record_name)
    rna_seq_record_name.match(SAMPLE_NAME_REGEX)[:sample_name]
  end
end
