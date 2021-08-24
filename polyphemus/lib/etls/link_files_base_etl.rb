require_relative "../metis_file_for_magma_model_etl"

class Polyphemus::LinkFilesBaseEtl < Polyphemus::MetisFileForMagmaModelEtl
  def initialize(attribute_name:, path_regex:, file_name_globs:, project_bucket_model_tuples:)
    @magma_models = {}
    @attribute_name = attribute_name
    super(
      project_bucket_model_tuples: project_bucket_model_tuples,
      file_name_globs: file_name_globs,
      metis_path_to_record_name_regex: path_regex,
    )
  end

  def process(cursor, files_by_record_name)
    logger.info("Processing files for: #{files_by_record_name.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: cursor[:project_name],
    )

    files_by_record_name.each do |file_record|
      next if should_skip_record?(file_record.record_name)
      next if file_record.files.empty?

      correct_record_name = corrected_record_name(link_model_record_name(file_record.files.first), file_record.record_name)
      update_request.update_revision(
        cursor[:model_name],
        correct_record_name,
        files_payload(cursor, file_record.files)
      )
      logger.info("Found #{file_record.files.length} files for #{correct_record_name}: #{file_record.files.map { |f| f.file_path }}")
    end

    magma_client.update_json(update_request)

    logger.info("Done")
  end

  private

  def magma_models(project_name)
    return @magma_models[project_name] if @magma_models.key?(project_name)

    @magma_models[project_name] = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
      project_name: project_name,
      model_name: "all",
      attribute_names: "all",
      record_names: [],
    )).models

    @magma_models[project_name]
  end

  def files_payload(cursor, files)
    {}.tap do |payload|
      payloads_for_attr = files.map do |file|
        serialize(cursor, file.file_path)
      end

      is_file_collection = is_file_collection?(cursor, @attribute_name)
      if payloads_for_attr.empty?
        payload[@attribute_name] = is_file_collection ? [] : blank_file
      else
        payload[@attribute_name] = is_file_collection ? payloads_for_attr : payloads_for_attr.first
      end
    end
  end

  def blank_file
    {
      "path": "::blank",
    }
  end

  def serialize(cursor, file_path)
    {
      path: metis_path(cursor, file_path),
      original_filename: ::File.basename(file_path),
    }
  end

  def metis_path(cursor, file_path)
    "metis://#{cursor[:project_name]}/#{cursor[:bucket_name]}/#{file_path}"
  end

  def is_file_collection?(cursor, attribute_name)
    magma_models(cursor[:project_name]).model(cursor[:model_name]).template.attributes.attribute(attribute_name.to_s).attribute_type == Etna::Clients::Magma::AttributeType::FILE_COLLECTION
  end

  def link_model_record_name(metis_file)
    raise "Should be implemented by subclasses"
  end

  def corrected_record_name(link_model_record_name, record_name)
    # Should be implemented by subclasses, if needed
    record_name
  end

  def should_skip_record?(record_name)
    # Should be implemented by subclasses, if needed
    false
  end
end
