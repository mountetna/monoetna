# Returns a set of Metis files, given a Magma model

require_relative "./metis_file_etl"

require_relative "etl"
require_relative "time_scan_based_etl_scanner"

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
    def initialize(project_bucket_pairs:, model_name: nil, metis_client: nil, magma_client: nil, limit: 20, file_name_globs: [])
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

      super(
        cursor_group: EtlCursorGroup.new(file_cursors),
        scanner: TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
          find_request = Etna::Clients::Metis::FindRequest.new(
            project_name: cursor[:project_name],
            bucket_name: cursor[:bucket_name],
          )

          prepare_base_find_request(cursor, find_request)

          find_request
        end.result_updated_at do |file_record|
          file_record.updated_at
        end.result_id do |file_record|
          file_record.file_path
        end.execute_batch_find do |find_request, i|
          records_request = magma_request(find_request.project_name)
          records_request.page_size = @limit * i
          # The limit will apply to the Magma client find.
          # Here we'll aggregate the Metis client results for each set
          #    of Magma client records. So in all likelihood
          #    we'll return more than @limit # of files.
          magma_record_names = self.magma_client.retrieve(
            records_request
          ).models.model(
            @model_name
          ).documents.document_keys

          all_files = []

          magma_record_names.each do |record_name|
            find_request_for_record = find_request.clone

            prepare_find_request_for_record(find_request_for_record, record_name)

            all_files.push(*self.metis_client.find(find_request_for_record).files.all)
          end

          all_files
        end,
      )
    end

    def magma_request(project_name)
      Etna::Clients::Magma::RetrievalRequest.new(
        project_name: project_name,
        model_name: @model_name,
        attribute_names: ["identifier"],
        record_names: "all",
        hide_templates: true,
        page: 1,
      )
    end

    def prepare_base_find_request(cursor, find_request)
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
        type: "file",
        attribute: "updated_at",
        predicate: ">=",
        value: (cursor.updated_at + 1).iso8601,
      )) unless cursor.updated_at.nil?
    end

    def prepare_find_request_for_record(find_request, record_name)
      # Only find files with the record_name in the path
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
        type: "file",
        attribute: "name",
        predicate: "glob",
        value: "#{record_name}/**/*",
      ))

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
  end
end
