require_relative 'etl'
require_relative 'etl_cursor'
require_relative 'time_scan_based_etl_scanner'

class Polyphemus
  class MagmaRecordEtlCursor < SimpleDelegator
    def initialize(job_name:, project_name:, model_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "model_name cannot be nil" if model_name.nil?
      super("#{job_name}_magma_records__#{project_name}_#{model_name}")
      self[:project_name] = project_name
      self[:model_name] = model_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  # Abstract base class for an ETL that scans Magma for records using the retrieve api.
  class MagmaRecordEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_model_pairs:, cursor_env:, scanner:, magma_client: nil, limit: 20, job_name: self.class.name, attribute_names: 'all')
      logger.info("Reading cursors...")

      cursors = cursors_from_pairs(
        pairs: project_model_pairs,
        pair_keys: %w[project_name model_name],
        cls: MagmaRecordEtlCursor,
        cursor_env: cursor_env
      )

      scanner = build_scanner if scanner.nil?

      @magma_client = magma_client
      @limit = limit
      @attribute_names = attribute_names

      super(cursors: cursors, scanner: scanner)
    end

    def build_scanner
      @scanner ||= TimeScanBasedEtlScanner.new.tap do |scanner|
        scanner.start_batch_state do |cursor|
          retrieve_request = Etna::Clients::Magma::RetrievalRequest.new(
            project_name: cursor[:project_name],
            model_name: cursor[:model_name],
            attribute_names: @attribute_names,
            order: 'updated_at',
            record_names: 'all',
            page: 1,
          )
          prepare_request(cursor, retrieve_request)
          retrieve_request
        end.result_updated_at do |record|
          updated_at = record.values.first['updated_at']
          updated_at ? Time.parse(updated_at) : Time.at(0)
        end.result_id do |record|
          record.keys.first
        end.execute_batch_find do |retrieve_request, i|
          retrieve_request.page_size = @limit * i
          documents = self.magma_client.retrieve(retrieve_request).models.model(retrieve_request.model_name).documents
          documents.document_keys.map { |i| { i => documents.document(i) } }
        end
      end
    end

    # Subclasses should override if they wish to adjust or add to the params of the retrieve request.
    def prepare_request(cursor, request)
      request.filter = "updated_at>=#{(cursor.updated_at + 1).iso8601}" unless cursor.updated_at.nil?

      if (end_at = cursor[:batch_end_at])
        request.filter = [
          request.filter,
          # This is crucial -- we use the <= here which can result in overlapping work executing because
          # we do not keep cursor state when running batches from the environment, and it is possible for a
          # process to run precisely at the time precision border between start and end times.
          # In this extreme edge case, prefer possibly duplication ensured completeness since jobs should be idempotent
          # anyways.
          "updated_at<=#{(end_at + 1).iso8601}"
        ]
      end
    end

    def find_batch
      super.map { |r| r.values.first }
    end
  end
end
