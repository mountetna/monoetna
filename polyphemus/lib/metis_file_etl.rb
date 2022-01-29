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
    def initialize(project_bucket_pairs:, metis_client: nil, limit: 20, file_name_params: {}, cursor_env: {}, scanner: build_scanner)
      @metis_client = metis_client
      @limit = limit
      @file_name_params = file_name_params

      cursors = cursors_from_pairs(
        pairs: project_bucket_pairs,
        pair_keys: %w[project_name bucket_name],
        cls: MetisFileEtlCursor,
        cursor_env: cursor_env
      )

      super(
        cursors: cursors,
        scanner: scanner,
      )
    end

    def serialize_batch(batch)
      super(batch.map(&:raw))
    end

    def deserialize_batch(string)
      super.map { |f| Etna::Clients::Metis::File.new(f) }
    end

    def build_scanner
      TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
        find_request = Etna::Clients::Metis::FindRequest.new(
          project_name: cursor[:project_name],
          bucket_name: cursor[:bucket_name],
          offset: 0
        )

        prepare_find_request(cursor, find_request)
        find_request
      end.result_updated_at do |file|
        file.updated_at
      end.result_id do |file|
        # This works as a change of the underlying file at a given name will trigger
        # further processing, and we're only concerned with the state of the nominal file
        # system at a point of time.
        file.file_path
      end.execute_batch_find do |find_request, i|
        execute_request(find_request, i)
      end
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

      if (end_at = cursor[:batch_end_at])
        find_request.add_param(Etna::Clients::Metis::FindParam.new(
          type: 'file',
          attribute: 'updated_at',
          predicate: '<=',
          value: (end_at + 1).iso8601,
        ))
      end

      begin
        @file_name_params.each do |predicate, values|
          values.each do |value|
            find_request.add_param(Etna::Clients::Metis::FindParam.new(
              type: 'file',
              attribute: 'name',
              predicate: predicate,
              value: value,
            ))
          end
        end
      end unless @file_name_params.keys.empty?
    end

    # Subclasses should override if they wish to filter or modify the output
    def execute_request(find_request, i)
      find_request.limit = @limit * i

      self.metis_client.find(find_request).files.all
    end
  end

  class MetisFileTailEtl < MetisFileEtl
  end
end
