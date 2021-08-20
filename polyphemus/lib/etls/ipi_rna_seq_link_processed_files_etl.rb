require_relative "../metis_file_for_magma_model_etl"
require_relative "../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqLinkProcessedFilesEtl < Polyphemus::MetisFileForMagmaModelEtl
  PATH_REGEX = /.*\/(?<record_name>.*)\/(?<original_file_name>.*\.(fastq\.gz|deduplicated\.cram|deduplicated\.cram\.crai|junction))$/
  PROJECT = "ipi"
  BUCKET = "data"
  MODEL = "rna_seq"
  CRAM_REGEX = /.*\.deduplicated\.cram$/
  CRAM_INDEX_REGEX = /.*\.deduplicated\.cram\.crai$/
  JUNCTION_REGEX = /.*\.junction$/
  UNMAPPED_FASTQ_REGEX = /.*\.fastq\.gz$/

  def initialize
    @helper = IpiHelper.new
    @magma_attributes = {}
    super(
      project_bucket_model_tuples: [[PROJECT, BUCKET, MODEL]],
      file_name_globs: ["bulkRNASeq/**/*", "output/**/*"],
      metis_path_to_record_name_regex: PATH_REGEX,
    )
  end

  def process(cursor, files_by_record_name)
    logger.info("Processing files for: #{files_by_record_name.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: cursor[:project_name],
    )

    files_by_record_name.each do |file_record|
      update_request.update_revision(
        cursor[:model_name],
        file_record.record_name,
        files_payload(cursor[:project_name], file_record.files)
      )
      logger.info("Found #{file_record.files.length} files for #{file_record.record_name}: #{file_record.files.map { |f| f.file_path }}")
    end

    magma_client.update_json(update_request)

    logger.info("Done")
  end

  def magma_attributes(project_name)
    return @magma_attributes[project_name] if @magma_attributes.key?(project_name)

    @magma_attributes[project_name] = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
      project_name: project_name,
      model_name: MODEL,
      attribute_names: attribute_regex.keys,
      record_names: [],
    )).models.model(MODEL).template.attributes

    @magma_attributes[project_name]
  end

  def files_payload(project_name, files)
    {}.tap do |payload|
      attribute_regex.each do |attribute_name, regex|
        file_payloads = files.select do |file|
          file.file_name =~ regex
        end.map do |file|
          serialize(file.file_path)
        end

        is_file_collection = is_file_collection?(project_name, attribute_name)
        if file_payloads.empty?
          payload[attribute_name] = is_file_collection ? [] : blank_file
        else
          payload[attribute_name] = is_file_collection ? file_payloads : file_payloads.first
        end
      end
    end
  end

  def blank_file
    {
      "path": "::blank",
    }
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

  def attribute_regex
    {
      "cram": CRAM_REGEX,
      "cram_index": CRAM_INDEX_REGEX,
      "junction": JUNCTION_REGEX,
      "unmapped_fastqs": UNMAPPED_FASTQ_REGEX,
    }
  end

  def is_file_collection?(project_name, attribute_name)
    magma_attributes(project_name).attribute(attribute_name.to_s).attribute_type == Etna::Clients::Magma::AttributeType::FILE_COLLECTION
  end
end
