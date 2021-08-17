require_relative "../metis_folder_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiCreateRnaSeqAndPlateRecordsEtl < Polyphemus::MetisFolderEtl
  PATH_REGEX = /.*\/(?<plate>plate\d+)_.*\/output\/(?<record_name>.*)/
  SAMPLE_NAME_REGEX = /^(?<sample_name>IPI.*\.[A-Z]+\d)\..*/
  NON_CANCER_REGEX = /(NASH|NAFLD)/
  PROJECT = "ipi"
  BUCKET = "data"

  def initialize
    @helper = IpiHelper.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_name_globs: ["output/*", "bulkRNASeq/*"],
    )
  end

  def process(cursor, folders)
    # All of these folders should be rna_seq record names.
    # We can extract the plate# from the folder_path:
    #   bulkRNASeq/processed/plate1_rnaseq_new/output/blahblahrecordname
    # We will also have to format the Control record names to match IPI validation.

    folders_by_plate = folders.group_by do |folder|
      plate(folder.folder_path)
    end

    logger.info("Ensuring plates exist: #{folders_by_plate.keys.join(",")}")
    ensure_plates(folders_by_plate.keys)

    create_records(folders_by_plate)

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

  def ensure_plates(plate_names)
    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )
    plate_names.each do |plate_name|
      update_request.update_revision("rna_seq_plate", plate_name, {
        "project" => "UCSF Immunoprofiler",
      })
    end
    magma_client.update_json(update_request)
  end

  def create_records(folders_by_plate)
    folders_by_plate.each do |plate_name, folders|
      update_request = Etna::Clients::Magma::UpdateRequest.new(
        project_name: PROJECT,
      )
      folders.each do |folder|
        next if is_non_cancer_sample?(folder.folder_name)

        record_name = @helper.is_control?(folder.folder_name) ?
          @helper.control_name(folder.folder_name) :
          @helper.corrected_rna_seq_tube_name(plate_name, folder.folder_name)

        attrs = {
          rna_seq_plate: plate_name,
        }

        attrs[:sample] = sample_name(record_name) unless @helper.is_control?(folder.folder_name)

        update_request.update_revision("rna_seq", record_name, attrs)
      end

      if update_request.revisions["rna_seq"]
        logger.info("Creating rna_seq records: #{update_request.revisions["rna_seq"].keys.join(",")}")
        magma_client.update_json(update_request)
      end
    end
  end

  def is_non_cancer_sample?(folder_name)
    # per Vincent
    # I guess the two major non-cancer samples that were gathered in IPIv1 were NAFLD/NASH and PSC (primary sclerosing cholangitis)
    folder_name =~ NON_CANCER_REGEX
  end

  def sample_name(rna_seq_record_name)
    rna_seq_record_name.match(SAMPLE_NAME_REGEX)[:sample_name]
  end
end
