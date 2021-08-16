require_relative "../metis_file_for_magma_model_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqLinkFastQsEtl < Polyphemus::MetisFileForMagmaModelEtl
  PATH_REGEX = /.*\/(?<record_name>.*)\/(?<original_file_name>.*\.fastq\.gz)$/
  PROJECT = "ipi"
  BUCKET = "integral_data"
  MODEL = "rna_seq"

  def initialize
    @helper = IpiHelper.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      file_name_globs: ["BulkRNASeq/**/*.fastq.gz"],
      model_name: MODEL,
      metis_path_to_record_name_regex: PATH_REGEX,
    )
  end

  def process(cursor, files_by_record_name)
    logger.info("Processing files for: #{files_by_record_name.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: PROJECT,
    )

    files_by_record_name.each do |file_record|
      update_request.update_revision("rna_seq", file_record.record_name, {
        "raw_fastq_files": files_payload(file_record.files),
      })
      logger.info("Found #{file_record.files.length} files for #{file_record.record_name}: #{file_record.files.map { |f| f.file_path }}")
    end

    magma_client.update_json(update_request)

    logger.info("Done")
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
