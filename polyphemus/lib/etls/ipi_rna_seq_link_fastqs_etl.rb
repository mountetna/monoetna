require_relative "../metis_file_for_magma_model_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqLinkFastQsEtl < Polyphemus::MetisFileForMagmaModelEtl
  PATH_REGEX = /.*\/(?<record_name>.*)\/(?<original_file_name>.*\.fastq\.gz)$/
  SAMPLE_NAME_REGEX = /^(?<sample_name>IPI.*\.[A-Z]+\d)\..*/
  PATIENT_IPI_NUMBER_REGEX = /^(?<ipi_number>IPI.*)\.[A-Z]+\d\..*/
  NON_CANCER_REGEX = /(NASH|NAFLD)/
  PROJECT = "ipi"
  BUCKET = "integral_data"
  MODEL = "rna_seq"

  def initialize
    @helper = IpiHelper.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      file_name_globs: ["BulkRNASeq/**/*.fastq.gz"],
      model_name: MODEL,
    )
  end

  def process(cursor, files)
    files_by_record_name = files.group_by do |file|
      record_name(file.file_path)
    end

    logger.info("Processing files: #{files.map { |f| f.file_path }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )

    files_by_record_name.each do |record_name, files|
      update_request.update_revision("rna_seq", record_name, {
        "raw_fastq_files": files_payload(files),
      })
    end

    magma_client.update_json(update_request)

    logger.info("Done")
  end

  def record_name(file_path)
    file_path.match(PATH_REGEX)[:record_name]
  end

  def files_payload(files)
    files.map do |file|
      serialize(file.file_path)
    end
  end

  def serialize(file_path)
    {
      path: metis_path(file_path),
      original_filename: ::File.basename(file_path),
    }
  end

  def metis_path(file_path)
    "metis://#{PROJECT}/#{BUCKET}/#{file_path}"
  end
end
