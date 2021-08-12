require_relative "../metis_folder_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqLinkFastQsEtl < Polyphemus::MetisFileEtl
  PATH_REGEX = /.*\/(?<record_name>.*)\/(?<original_file_name>.*\.fastq\.gz)$/
  SAMPLE_NAME_REGEX = /^(?<sample_name>IPI.*\.[A-Z]+\d)\..*/
  PATIENT_IPI_NUMBER_REGEX = /^(?<ipi_number>IPI.*)\.[A-Z]+\d\..*/
  NON_CANCER_REGEX = /(NASH|NAFLD)/
  PROJECT = "ipi"
  BUCKET = "integral_data"

  def initialize
    @helper = IpiHelper.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      file_name_globs: ["BulkRNASeq/**/*.fastq.gz"],
    )
  end

  def process(cursor, files)
    # We'll get all the raw fastq files, and then match them up with Magma records.
    #   If there is no corresponding Magma record, do nothing.
    # We can't do an "OR" find in Metis yet, so we match up to records this way,
    #   instead of scanning for rna_seq records and then building Metis queries.

    files_by_record_name = files.group_by do |file|
      record_name(file.file_path)
    end

    logger.info("Processing files: #{files.map { |f| f.file_path }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )
    existing_magma_records = magma_records(files_by_record_name.keys)
    existing_magma_records.document_keys.each do |record_name|
      record = existing_magma_records.document(record_name)
      raw_fastq_files = record["raw_fastq_files"].nil? ? [] : record["raw_fastq_files"]

      files_by_record_name[record_name].each do |new_fastq_file|
        raw_fastq_files << serialize(new_fastq_file.file_path)
      end

      update_request.update_revision("rna_seq", record_name, {
        "raw_fastq_files": raw_fastq_files,
      })
    end

    magma_client.update_json(update_request)

    logger.info("Done")
  end

  def record_name(file_path)
    file_path.match(PATH_REGEX)[:record_name]
  end

  def magma_records(potential_record_names)
    # This should always be a subset of files_by_record_name.keys, and is
    #   a check that the Magma record exists.
    magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
      project_name: PROJECT,
      model_name: "rna_seq",
      attribute_names: ["tube_name", "raw_fastq_files"],
      record_names: potential_record_names,
    )).models.model("rna_seq").documents
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
