require_relative "../helpers"

class Polyphemus::MetisFilesLinkerBase
  include WithLogger
  include WithEtnaClients

  attr_reader :bucket_name, :project_name

  def initialize(project_name:, bucket_name:)
    @project_name = project_name
    @bucket_name = bucket_name
    @magma_models = {}
  end

  # Important!  If files_by_record_name includes files to be included in a file collection, /all files belonging to that collection must be present/.
  # the metis file in watch folder etl uses folder_ids_matching_file_collections to grab all containing folders for files that
  # match file collections and ensures that a folder list is performed and that all known files in that grouping are included.
  def link(model_name:, files_by_record_name:, attribute_regex:)
    return false if files_by_record_name.empty?

    logger.info("Processing files for: #{files_by_record_name.map { |f| f.record_name }.join(",")}")

    update_request = Etna::Clients::Magma::UpdateRequest.new(
      project_name: project_name,
    )

    # Currently, there isn't a good magma api for transactions or deltas applied to collections, so this
    # approach is vulnerable to race conditions such as a race to override by two simult writers.
    # In this case, we should be ok, given the assumption that we only run one instance of the ETL or shard by record id,
    # and that we don't plan on mixing source of truth (linking and hand appending files).
    record_batch_ids = files_by_record_name.map(&:record_name)
    record_batch = magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
      project_name: project_name,
      model_name: model_name,
      page_size: record_batch_ids.length, # Grab everything in one go, this should still have a fairly low batched ceiling,
      record_names: record_batch_ids,
      hide_templates: true,
    )).models.model(model_name)&.documents

    if record_batch.nil?
      raise "Unexpect nil in retrieval documents response for linker.  Did the server response change?"
    end

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
        record: record_batch.document(file_record.record_name),
      )

      update_request.update_revision(
        model_name,
        correct_record_name,
        payload
      )
      logger.info("Found #{payload.values.flatten.length} files for #{correct_record_name}: #{payload.values.flatten.join(",")}")
    end

    magma_client.update_json(update_request) if update_request.revisions.keys.length > 0

    logger.info("Done")

    update_request.revisions
  end

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

  def files_payload(project_name:, model_name:, files:, attribute_regex:, record:)
    {}.tap do |payload|
      attribute_regex.each do |attribute_name, regex|
        payloads_for_attr = files.select do |file|
          file.file_name =~ regex
        end.map do |file|
          serialize(file)
        end

        next if payloads_for_attr.empty?

        if is_file_collection?(project_name, model_name, attribute_name)
          if record.nil? || (existing_files = record[attribute_name]).nil?
            payload[attribute_name] = payloads_for_attr
          else
            payload_by_paths = payloads_for_attr.map { |f| f[:path] }
            payload[attribute_name] = payloads_for_attr + existing_files.select do |file|
              !payload_by_paths.include?(file[:path])
            end
          end
        else
          payload[attribute_name] = payloads_for_attr.first
        end
      end
    end
  end

  def serialize(file)
    {
      path: metis_path(file),
      original_filename: file.file_name,
    }
  end

  def current_magma_record_names(project_name, model_name)
    request = Etna::Clients::Magma::RetrievalRequest.new(
      project_name: project_name,
      model_name: model_name,
      attribute_names: ["identifier"],
      record_names: "all",
      hide_templates: true,
    )
    self.magma_client.retrieve(request).models.model(model_name).documents.document_keys
  end

  def metis_path(file)
    # Technically you could build an ETL to watch for a file without a folder...
    #   but that doesn't seem like a realistic use case.
    "metis://#{file.project_name}/#{file.bucket_name}/#{watch_folder_for_file(file).folder_path}/#{file.file_name}"
  end

  def watch_folder_for_file(file)
    Polyphemus::WatchFolder.where(
      project_name: file.project_name,
      bucket_name: file.bucket_name,
      watch_type: "link_files",
      metis_id: file.folder_id,
    ).first
  end

  def full_path_for_file(file)
    watch_folder = watch_folder_for_file(file)

    unless watch_folder.nil?
      "#{watch_folder.folder_path}/#{file.file_name}"
    end
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

  def organize_metis_files_by_magma_record(
    metis_files:,
    magma_record_names:,
    path_regex:,
    record_name_gsub_pair: nil
  )
    metis_files_by_record_name = metis_files.group_by do |file|
      match = full_path_for_file(file)&.match(path_regex)

      if match
        record_name = corrected_record_name(match[:record_name])

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
