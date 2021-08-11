require_relative "../metis_folder_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiCreateRnaSeqAndPlateRecordsEtl < Polyphemus::MetisFolderEtl
  PATH_REGEX = /.*\/(?<plate>plate\d+)_.*\/output\/(?<record_name>.*)/
  SAMPLE_NAME_REGEX = /^(?<sample_name>IPI.*\.[A-Z]+\d)\..*/
  PATIENT_IPI_NUMBER_REGEX = /^(?<ipi_number>IPI.*)\.[A-Z]+\d\..*/
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

    logger.error("Unmatched folder: #{folder_path}") unless match

    match[:plate].capitalize
  end

  def ensure_plates(plate_names)
    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )
    plate_names.each do |plate_name|
      update_request.update_revision("rna_seq_plate", plate_name, {})
    end
    magma_client.update_json(update_request)
  end

  def create_records(folders_by_plate)
    folders_by_plate.each do |plate_name, folders|
      update_request = Etna::Clients::Magma::UpdateRequest.new(
        project_name: PROJECT,
      )
      folders.each do |folder|
        record_name = is_control?(folder.folder_name) ?
          control_name(folder.folder_name) :
          folder.folder_name

        attrs = {
          rna_seq_plate: plate_name,
        }

        if !is_control?(folder.folder_name)
          attrs[:sample] = sample_name(record_name)
          containing_record_workflow.ensure_record("sample", containing_records(record_name))
        end

        update_request.update_revision("rna_seq", record_name, attrs)
      end
      logger.info("Creating rna_seq records: #{update_request.revisions["rna_seq"].keys.join(",")}")
      magma_client.update_json(update_request)
    end
  end

  def is_control?(record_name)
    record_name =~ /^control/i
  end

  def control_name(folder_name)
    # Control_(UHR|Jurkat).Plate\d+
    control, plate = folder_name.split(".")
    _, control_type = control.split("_")

    "Control_#{control_type =~ /jurkat/i ? "Jurkat" : "UHR"}.#{plate.capitalize}"
  end

  def remove_extraneous_sample_letters(sample_name)
    # Because of IPIHEP033.NASH1 -> IPIHEP033.N1 on plate 26
    parts = sample_name.split(".")
    [parts.first, parts.last[0] + parts.last[-1]].join(".")
  end

  def sample_name(rna_seq_record_name)
    remove_extraneous_sample_letters(
      rna_seq_record_name.match(SAMPLE_NAME_REGEX)[:sample_name]
    )
  end

  def patient_ipi_number(rna_seq_record_name)
    rna_seq_record_name.match(PATIENT_IPI_NUMBER_REGEX)[:ipi_number]
  end

  def magma_crud
    @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
      magma_client: magma_client,
      project_name: PROJECT,
      read_only: false,
    )
  end

  def containing_records(record_name)
    {
      "sample" => sample_name(record_name),
      "patient" => patient_ipi_number(record_name),
      "experiment" => @helper.experiment_from_patient_number(patient_ipi_number(record_name)),
      "project" => "UCSF Immunoprofiler",
    }
  end

  def containing_record_workflow
    @containing_record_workflow ||= Etna::Clients::Magma::EnsureContainingRecordWorkflow.new(magma_crud: magma_crud, models: models)
  end

  def models
    @models ||= begin
        magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(project_name: PROJECT, model_name: "all")).models
      end
  end
end
