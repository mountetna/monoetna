require_relative 'etl'
require_relative 'time_scan_based_etl_scanner'

class Polyphemus
  class MagmaRecordEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, model_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "model_name cannot be nil" if model_name.nil?
      super("#{job_name}_magma_records__#{project_name}_#{model_name}")
      self[:project_name] = project_name
      self[:model_name] = model_name
    end

    def reset!
      self[:seen_ids] = []
      super
    end
  end

  # Abstract base class for an ETL that scans metis for files using the find api.
  class MagmaFileEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_model_pairs:, magma_client: nil, limit: 20)
      cursors = project_model_pairs.map do |project_name, model_name|
        MagmaRecordEtlCursor.new(job_name: self.class.name, project_name: project_name, model_name: model_name).load_from_db
      end

      @magma_client = magma_client
      @limit = limit

      super(
          cursor_group: EtlCursorGroup.new(cursors),
          scanner: TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
            find_request = Etna::Clients::Metis::FindRequest.new(
                project_name: cursor[:project_name],
                bucket_name: cursor[:bucket_name],
            )
            prepare_find_request(cursor, find_request)
            find_request
          end.result_updated_at do |file|
            file.updated_at
          end.result_id do |file|
            file.file_path
          end.execute_batch_find do |find_request, offset|
            find_request.offset = offset
            metis_client.find(find_request).files.all
          end
      )
    end

    # Subclasses should override if they wish to adjust or add to the params of the find request.
    def prepare_find_request(cursor, find_request)
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
              type: 'file',
              attribute: 'updated_at',
              predicate: '>=',
              value: (cursor.updated_at + 1).iso8601,
          )
      ) unless cursor.updated_at.nil?
      find_request.limit = @limit
    end
  end
end
