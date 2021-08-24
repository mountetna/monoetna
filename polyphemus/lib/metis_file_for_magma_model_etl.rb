# Returns a set of Metis files, given a Magma model

require_relative "./metis_file_etl"

require_relative "etl"
require_relative "hash_scan_based_etl_scanner"

class Polyphemus
  class MetisFileForMagmaModelEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, bucket_name:, model_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      raise "model_name cannot be nil" if model_name.nil?
      super("#{job_name}_metis_files_for_magma_model_#{project_name}_#{bucket_name}_#{model_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
      self[:model_name] = model_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  class MetisFileForMagmaModelEtl < Etl
    def initialize(
      project_bucket_model_tuples:,
      metis_client: nil,
      magma_client: nil,
      file_name_globs: [],
      metis_path_to_record_name_regex: nil,
      record_name_gsub_pair: nil,
      limit: 20
    )
      file_cursors = project_bucket_model_tuples.map do |project_name, bucket_name, model_name|
        MetisFileForMagmaModelEtlCursor.new(
          job_name: self.class.name,
          project_name: project_name,
          bucket_name: bucket_name,
          model_name: model_name,
        ).load_from_db
      end

      @metis_client = metis_client
      @magma_client = magma_client
      @file_name_globs = file_name_globs
      @path_regex = metis_path_to_record_name_regex
      @record_name_gsub_pair = record_name_gsub_pair
      @limit = limit

      @metis_file_limit = 200
      @metis_files_for_cursor = {}

      super(
        cursor_group: EtlCursorGroup.new(file_cursors),
        scanner: HashScanBasedEtlScanner.new.start_batch_state do |cursor|
          calculate_metis_files_for_cursor(cursor)
        end.result_updated_at do |metis_files_record|
          metis_files_record.updated_at
        end.result_file_hashes do |metis_files_record|
          # Only record a changed set if the hashes or files change
          #   for a given Magma record
          metis_files_record.file_paths_hashes
        end.result_id do |metis_files_record|
          metis_files_record.record_name
        end.execute_batch_find do |metis_files_for_magma_records, i|
          metis_files_for_magma_records.slice(@limit * (i - 1), @limit) || []
        end,
      )
    end

    def calculate_metis_files_for_cursor(cursor)
      cursor_string_sym = "#{cursor[:project_name]}_#{cursor[:bucket_name]}_#{cursor[:model_name]}".to_sym

      return @metis_files_for_cursor[cursor_string_sym] if @metis_files_for_cursor.has_key?(cursor_string_sym)

      metis_request = Etna::Clients::Metis::FindRequest.new(
        project_name: cursor[:project_name],
        bucket_name: cursor[:bucket_name],
        offset: 0,
        limit: @metis_file_limit,
      )

      magma_request = Etna::Clients::Magma::RetrievalRequest.new(
        project_name: cursor[:project_name],
        model_name: cursor[:model_name],
        attribute_names: ["identifier"],
        record_names: "all",
        hide_templates: true,
      )

      magma_record_names = fetch_magma_record_names(magma_request)

      prepare_metis_find_request(metis_request)

      files = magma_record_names.empty? ?
        [] :
        organize_metis_files_by_magma_record(
        collect_all_metis_files(metis_request),
        magma_record_names
      )

      @metis_files_for_cursor[cursor_string_sym] = files

      files
    end

    def prepare_metis_find_request(find_request)
      # Do not track the cursor's updated_at, because
      #   we always need to get all files on Metis that match the given globs,
      #   so we can link all files to Magma records, i.e. for a file_collection,
      #   if a file is added later, we need it to be found with all the previous
      #   files.
      begin
        @file_name_globs.each do |glob|
          find_request.add_param(Etna::Clients::Metis::FindParam.new(
            type: "file",
            attribute: "name",
            predicate: "glob",
            value: glob,
          ))
        end
      end unless @file_name_globs.empty?
    end

    def fetch_magma_record_names(magma_request)
      self.magma_client.retrieve(
        magma_request
      ).models.model(
        magma_request.model_name
      ).documents.document_keys
    end

    def collect_all_metis_files(metis_request)
      files = []
      i = 0

      loop do
        metis_request.offset = @metis_file_limit * i
        new_files = self.metis_client.find(metis_request).files.all

        break if new_files.empty?

        files.push(*new_files)
        i += 1
      end

      files
    end

    def organize_metis_files_by_magma_record(metis_files, magma_record_names)
      metis_files_by_record_name = metis_files.group_by do |file|
        match = file.file_path.match(@path_regex)

        if match
          record_name = match[:record_name]

          record_name = record_name.gsub(@record_name_gsub_pair.first, @record_name_gsub_pair.last) if @record_name_gsub_pair

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
