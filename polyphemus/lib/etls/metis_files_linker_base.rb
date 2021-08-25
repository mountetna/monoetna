class Polyphemus::MetisFilesLinkerBase
  attr_reader :logger

  def initialize(logger)
    @magma_models = {}
    @logger = logger
  end

  def process(project_name:, model_name:, files_by_record_name:, attribute_regex:)
    logger.info("Processing files for: #{files_by_record_name.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: project_name,
    )

    files_by_record_name.each do |file_record|
      next if should_skip_record?(file_record.record_name)
      next if file_record.files.empty?

      correct_record_name = corrected_record_name(
        link_model_record_name(file_record.files.first),
        file_record.record_name
      )
      update_request.update_revision(
        model_name,
        correct_record_name,
        files_payload(project_name, model_name, file_record.files)
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

  def files_payload(project_name, model_name, files)
    {}.tap do |payload|
      payloads_for_attr = files.map do |file|
        serialize(file)
      end

      is_file_collection = is_file_collection?(project_name, model_name, @attribute_name)

      payload[@attribute_name] = is_file_collection ? payloads_for_attr : payloads_for_attr.first
    end
  end

  def serialize(file)
    {
      path: metis_path(file),
      original_filename: ::File.basename(file.file_path),
    }
  end

  def metis_path(file)
    "metis://#{file.project_name}/#{file.bucket_name}/#{file.file_path}"
  end

  def is_file_collection?(project_name, model_name, attribute_name)
    magma_models(project_name).model(model_name).template.attributes.attribute(attribute_name.to_s).attribute_type == Etna::Clients::Magma::AttributeType::FILE_COLLECTION
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
