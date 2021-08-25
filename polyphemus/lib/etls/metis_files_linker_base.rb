require_relative "../helpers"

class Polyphemus::MetisFilesLinkerBase
  include WithLogger
  include WithEtnaClients

  def initialize
    @magma_models = {}
  end

  def link(project_name:, model_name:, files_by_record_name:, attribute_regex:)
    logger.info("Processing files for: #{files_by_record_name.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: project_name,
    )

    files_by_record_name.each do |file_record|
      next if should_skip_record?(file_record.record_name)
      next if file_record.files.empty?

      correct_record_name = corrected_record_name(
        file_record.record_name
      )

      payload = files_payload(
        project_name: project_name,
        model_name: model_name,
        files: file_record.files,
        attribute_regex: attribute_regex,
      )

      update_request.update_revision(
        model_name,
        correct_record_name,
        payload
      )
      logger.info("Found #{payload.values.flatten.length} files for #{correct_record_name}: #{payload.values.flatten.join(",")}")
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

  def files_payload(project_name:, model_name:, files:, attribute_regex:)
    {}.tap do |payload|
      attribute_regex.each do |attribute_name, regex|
        payloads_for_attr = files.select do |file|
          file.file_name =~ regex
        end.map do |file|
          serialize(file)
        end

        next if payloads_for_attr.empty?

        payload[attribute_name] = is_file_collection?(project_name, model_name, attribute_name) ?
          payloads_for_attr :
          payloads_for_attr.first
      end
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

  def corrected_record_name(record_name)
    # Should be implemented by subclasses, if needed
    record_name
  end

  def should_skip_record?(record_name)
    # Should be implemented by subclasses, if needed
    false
  end

  def self.organize_metis_files_by_magma_record(
    metis_files:,
    magma_record_names:,
    path_regex:,
    record_name_gsub_pair: nil
  )
    metis_files_by_record_name = metis_files.group_by do |file|
      match = file.file_path.match(path_regex)

      if match
        record_name = match[:record_name]

        record_name = record_name.gsub(record_name_gsub_pair.first, record_name_gsub_pair.last) if record_name_gsub_pair

        record_name
      else
        nil
      end
    end

    magma_record_names.map do |magma_record_name|
      next if metis_files_by_record_name[magma_record_name].nil?

      MetisFilesForMagmaRecord.new(
        magma_record_name,
        metis_files_by_record_name[magma_record_name]
      )
    end.compact
  end

  class MetisFilesForMagmaRecord
    attr_reader :record_name, :files

    def initialize(magma_record_name, metis_files)
      @record_name = magma_record_name
      @files = metis_files
    end

    def file_paths_hashes
      files.map do |file|
        [file.file_path, file.file_hash]
      end
    end

    def updated_at
      files.map do |file|
        file.updated_at
      end.minmax.last
    end
  end
end
