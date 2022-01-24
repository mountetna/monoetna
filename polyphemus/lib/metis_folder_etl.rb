require_relative "etl"
require_relative "time_scan_based_etl_scanner"

class Polyphemus
  class MetisFolderEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      super("#{job_name}_metis_folders_#{project_name}_#{bucket_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  # Abstract base class for an ETL that scans metis for folders using the find api.
  class MetisFolderEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, metis_client: nil, cursor_env:, scanner:, folder_name_regexes: [], folder_name_globs: [], limit: 20)
      @metis_client = metis_client
      @limit = limit
      @folder_name_globs = folder_name_globs
      @folder_name_regexes = folder_name_regexes

      cursors = cursors_from_pairs(
        pairs: project_bucket_pairs,
        pair_keys: %w[project_name bucket_name],
        cls: MetisFolderEtlCursor,
        cursor_env: cursor_env
      )

      scanner = build_scanner if scanner.nil?

      super(
        cursors: cursors,
        scanner: scanner
      )
    end

    def build_scanner
      TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
        find_request = Etna::Clients::Metis::FindRequest.new(
          project_name: cursor[:project_name],
          bucket_name: cursor[:bucket_name],
        )
        prepare_find_request(cursor, find_request)
        find_request
      end.result_updated_at do |folder|
        folder.updated_at
      end.result_id do |folder|
        folder.folder_path
      end.execute_batch_find do |find_request, i|
        find_request.limit = @limit * i
        self.metis_client.find(find_request).folders.all
      end
    end

    # Subclasses should override if they wish to adjust or add to the params of the find request.
    def prepare_find_request(cursor, find_request)
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
        type: "folder",
        attribute: "updated_at",
        predicate: ">=",
        value: (cursor.updated_at + 1).iso8601,
      )) unless cursor.updated_at.nil?

      if (end_at = cursor[:batch_end_at])
        find_request.add_param(Etna::Clients::Metis::FindParam.new(
          type: 'folder',
          attribute: 'updated_at',
          predicate: '<=',
          value: (end_at + 1).iso8601,
        ))
      end

      begin
        @folder_name_globs.each do |glob|
          find_request.add_param(Etna::Clients::Metis::FindParam.new(
            type: "folder",
            attribute: "name",
            predicate: "glob",
            value: glob,
          ))
        end
      end unless @folder_name_globs.empty?
      begin
        @folder_name_regexes.each do |regex|
          find_request.add_param(Etna::Clients::Metis::FindParam.new(
            type: "folder",
            attribute: "name",
            predicate: "=~",
            value: regex,
          ))
        end
      end unless @folder_name_regexes.empty?
    end
  end
end
