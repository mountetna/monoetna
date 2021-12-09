require_relative "../helpers"

class Polyphemus::MetisFilesLinkerBase
  include WithLogger
  include WithEtnaClients

  attr_reader :bucket_name, :project_name

  DOWNLOAD_REGEX = /^https:\/\/.*\/(?<project_name>.*)\/download\/(?<bucket_name>.*)\/(?<file_path>[^\?]*).*$/

  def initialize(project_name:, bucket_name:, record_name_regex:)
    @project_name = project_name
    @bucket_name = bucket_name
    @magma_models = {}
    @record_name_regex = record_name_regex
  end

  def cursor_model_group_state(cursor, model_name, group_name, default)
    linked_models_state = cursor[:linked_models_state] ||= {}
    model_state = linked_models_state[model_name.to_s] ||= {}
    model_state[group_name.to_s] ||= default
  end

  # Important!  If files_by_record_name includes files to be included in a file collection, /all files belonging to that collection must be present/.
  # the metis file in watch folder etl uses folder_ids_matching_file_collections to grab all containing folders for files that
  # match file collections and ensures that a folder list is performed and that all known files in that grouping are included.
  def link(model_name:, files_by_record_name:, attribute_regex:, cursor: nil)
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
        record_name: file_record.record_name,
        cursor: cursor
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

  def files_payload(project_name:, model_name:, files:, attribute_regex:, record:, cursor: nil, record_name:)
    {}.tap do |payload|
      attribute_regex.each do |attribute_name, regex|
        payloads_for_attr = files.select do |file|
          file.file_name =~ regex
        end.map do |file|
          serialize(file)
        end

        next if payloads_for_attr.empty?

        if is_file_collection?(project_name, model_name, attribute_name)
          if record.nil? || (existing_files = record[attribute_name.to_s]).nil?
            payload[attribute_name] = payloads_for_attr
          else
            payload_by_filename = payloads_for_attr.map { |f| f["original_filename"] }
            payload[attribute_name] = payloads_for_attr + convert_urls_to_metis_path(existing_files.select do |file|
              !payload_by_filename.include?(file["original_filename"])
            end)
          end
        else
          payload[attribute_name] = payloads_for_attr.first
        end
      end
    end
  end

  def serialize(file)
    {
      "path" => metis_path(file),
      "original_filename" => file.file_name,
    }
  end

  def convert_urls_to_metis_path(existing_files)
    #  Convert the download "url" value for existing files
    #   to a metis://<path> so that we can re-send it back
    #   to Magma in an update. Otherwise Magma doesn't know
    #   what to do with a download URL and non-metis path.
    existing_files.map do |file|
      {
        "original_filename" => file["original_filename"],
        "path" => metis_path_from_download_url(file["url"])
      }
    end
  end

  def metis_path_from_download_url(download_url)
    match = download_url.match(DOWNLOAD_REGEX)

    begin
      logger.error("Found an invalid existing file #{download_url}")
      return nil
    end if match.nil?

    "metis://#{match[:project_name]}/#{match[:bucket_name]}/#{match[:file_path]}"
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
    "metis://#{file.project_name}/#{file.bucket_name}/#{file.file_path}"
  end

  def is_file_collection?(project_name, model_name, attribute_name)
    attribute = magma_models(project_name).model(model_name).template.attributes.attribute(attribute_name.to_s)
    raise "Unknown linker attribute #{attribute_name}" if attribute.nil?
    attribute.attribute_type == Etna::Clients::Magma::AttributeType::FILE_COLLECTION
  end

  def corrected_record_name(record_name)
    # Should be implemented by subclasses, if needed
    record_name
  end

  def should_skip_record?(record_name)
    # Should be implemented by subclasses, if needed
    false
  end

  def record_name_by_path(file_or_folder)
    if file_or_folder.is_a?(Etna::Clients::Metis::File)
      match = file_or_folder.file_path.match(@record_name_regex)
    elsif file_or_folder.is_a?(Etna::Clients::Metis::Folder)
      match = (file_or_folder.folder_path + '/').match(@record_name_regex)
    else
      raise "file_or_folder must be either a File or Folder object"
    end

    if match
      corrected_record_name(match[:record_name])
    else
      nil
    end
  end

  def organize_metis_files_by_magma_record(
    metis_files:,
    magma_record_names:
  )
    metis_files_by_record_name = metis_files.group_by do |file|
      next if file.file_path.nil?
      record_name_by_path(file)
    end

    metis_files_by_record_name.keys.map do |matched_record_name|
      unless magma_record_names.include?(matched_record_name)
        logger.info("Could not find magma record matching file record #{matched_record_name}")
        next
      end

      MetisFilesForMagmaRecord.new(
        matched_record_name,
        metis_files_by_record_name[matched_record_name]
      )
    end.compact
  end

  class MetisFilesForMagmaRecord
    attr_reader :record_name, :files

    def initialize(record_name, metis_files)
      @record_name = record_name
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
