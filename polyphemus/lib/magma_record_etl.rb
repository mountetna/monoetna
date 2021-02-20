require_relative 'etl'
require_relative 'etl_cursor'
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
      super { self[:seen_ids] = [] }
    end
  end

  # Abstract base class for an ETL that scans metis for files using the find api.
  class MagmaRecordEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_model_pairs:, magma_client: nil, limit: 20, job_name: self.class.name)
      logger.info("Reading cursors...")
      cursors = project_model_pairs.map do |project_name, model_name|
        MagmaRecordEtlCursor.new(job_name: job_name, project_name: project_name, model_name: model_name).load_from_db
      end

      @magma_client = magma_client
      @limit = limit

      super(
          cursor_group: EtlCursorGroup.new(cursors),
          scanner: TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
            retrieve_request = Etna::Clients::Magma::RetrievalRequest.new(
                project_name: cursor[:project_name],
                model_name: cursor[:model_name],
                order: 'updated_at',
                record_names: 'all',
                page: 1,
            )
            prepare_retrieve_request(cursor, retrieve_request)
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
      )
    end

    # Subclasses should override if they wish to adjust or add to the params of the retrieve request.
    def prepare_retrieve_request(cursor, retrieve_request)
      retrieve_request.filter = "updated_at>=#{(cursor.updated_at + 1).iso8601}" unless cursor.updated_at.nil?
    end

    def find_batch
      super.map { |r| r.values.first }
    end
  end
end
