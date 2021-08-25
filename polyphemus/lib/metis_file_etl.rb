require_relative 'etl'
require_relative 'time_scan_based_etl_scanner'

class Polyphemus
  class MetisFileEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      super("#{job_name}_metis_files_#{project_name}_#{bucket_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  # Abstract base class for an ETL that scans metis for files using the find api.
  class MetisFileEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, metis_client: nil, limit: 20)
      file_cursors = project_bucket_pairs.map do |project_name, bucket_name|
        MetisFileEtlCursor.new(job_name: self.class.name, project_name: project_name, bucket_name: bucket_name).load_from_db
      end

      @metis_client = metis_client
      @limit = limit

      super(
          cursor_group: EtlCursorGroup.new(file_cursors),
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
          end.execute_batch_find do |find_request, i|
            find_request.limit = @limit * i
            execute_request(find_request)
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
    end

    # Subclasses should override if they wish to filter or modify the output
    def execute_request(find_request)
      self.metis_client.find(find_request).files.all
    end
  end
end
