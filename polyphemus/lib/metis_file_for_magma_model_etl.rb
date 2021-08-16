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
      project_bucket_pairs:,
      model_name: nil,
      metis_client: nil,
      magma_client: nil,
      limit: 20,
      file_name_globs: [],
      metis_path_to_record_name_regex: nil
    )
      file_cursors = project_bucket_pairs.map do |project_name, bucket_name|
        MetisFileForMagmaModelEtlCursor.new(
          job_name: self.class.name,
          project_name: project_name,
          bucket_name: bucket_name,
          model_name: model_name,
        ).load_from_db
      end

      @metis_client = metis_client
      @magma_client = magma_client
      @limit = limit
      @file_name_globs = file_name_globs
      @model_name = model_name
      @path_regex = metis_path_to_record_name_regex

      super(
        cursor_group: EtlCursorGroup.new(file_cursors),
        scanner: HashScanBasedEtlScanner.new.start_batch_state do |cursor|
          metis_request = Etna::Clients::Metis::FindRequest.new(
            project_name: cursor[:project_name],
            bucket_name: cursor[:bucket_name],
          )

          magma_request = Etna::Clients::Magma::RetrievalRequest.new(
            project_name: cursor[:project_name],
            model_name: @model_name,
            attribute_names: ["identifier"],
            record_names: "all",
            hide_templates: true,
          )

          record_names = fetch_magma_record_names(magma_request)

          prepare_metis_find_request(metis_request)

          # Pass both elements back instead of trying
          #   to include the record names in the Metis
          #   request directly, so that in
          #   execute_batch_find we can better handle
          #   cases when there are no Magma records.
          [metis_request, record_names]
        end.result_updated_at do |metis_files_record|
          metis_files_record.updated_at
        end.result_file_hashes do |metis_files_record|
          metis_files_record.file_paths_hashes
        end.result_id do |metis_files_record|
          metis_files_record.record_name
        end.execute_batch_find do |request_records_pair, i|
          metis_request = request_records_pair.first
          magma_record_names = request_records_pair.last

          return [] unless magma_record_names.length > 0

          metis_files = self.metis_client.find(metis_request).files.all
          metis_files_by_record_name = metis_files.group_by do |file|
            file.file_path.match(@path_regex)[:record_name]
          end

          magma_record_names.map do |magma_record_name|
            next if metis_files_by_record_name[magma_record_name].nil?

            MetisFilesForMagmaRecord.new(
              magma_record_name,
              metis_files_by_record_name[magma_record_name]
            )
          end.compact
        end,
      )
    end

    def prepare_metis_find_request(find_request)
      # Do not track the cursor's updated_at, because
      #   in order to find files that were deleted, we always need
      #   to get all files on Metis that match the given globs.
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
        @model_name
      ).documents.document_keys
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
