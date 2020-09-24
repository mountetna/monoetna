require 'ostruct'
require_relative 'helpers'

class Polyphemus
  class MetisFileEtlCursor < EtlCursor
    def initialize(project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      super("metis_file_etl_#{project_name}_#{bucket_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
    end

    def include?(file)
      return true if super(file.updated_at)
      ids = value[:ids]
      return false if ids.nil?
      ids.include?(file.id)
    end
  end

  # Abstract base class meant to do things.
  class MetisFileEtl < Etna::Command
    # Will use default client based on environment unless override is set in constructor
    include WithEtnaClients
    include WithLogger

    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, metis_client: nil)
      @metis_client = metis_client
      @cursor_group = EtlCursorGroup.new(project_bucket_pairs.map { |project_name, bucket_name| MetisFileEtlCursor.new(project_name: project_name, bucket_name: bucket_name).load_from_db })
    end

    # Subclasses should override if they wish to adjust or add to the params of the find request.
    def prepare_find_request(find_request)
      find_request.add_param(Etna::Clients::Metis::FindParam.new(
              type: 'file',
              attribute: 'updated_at',
              predicate: '>=',
              value: cursor.updated_at,
          )
      ) unless cursor.updated_at.nil?

      if cursor[:offset].nil? || cursor[:limit].nil?
        find_request.offset = 0
      else
        find_request.offset = cursor[:offset] + cursor[:limit]
      end

      find_request.limit = cursor[:limit] || 20
    end

    # Returns true iff some files were processed as part of this iteration.
    def run_once
      @cursor_group.with_next do |cursor|
        find_request = Etna::Clients::Metis::FindRequest.new(
            project_name: cursor[:project_name],
            bucket_name: cursor[:bucket_name],
        )
        prepare_find_request(find_request)

        response = metis_client.find(find_request)
        files = response.files.all.filter { |f| !cursor.include?(f) }
        return false if files.empty?

        process(cursor, files)
        update_cursor(cursor, find_request, files)
        return true
      end

      false
    end

    # Subclasses should override with their processing loop.
    # IMPORTANT!  Keep in mind that because the process can fail part way through work, there is no guarantee of single
    # delivery and that process should be idempotent to the best of its ability.
    # This class does best effort to manage work 'split' along the same updated at to try to ensure that for any
    # successful process and update_cursor call, a duplicate one should not occur after.  However this is only an
    # optimization and not a guarantee.
    def process(cursor, files)
    end

    def update_cursor(cursor, file_request, processed_files)
      last_file = processed_files.max { |a, b| [a.updated_at, a.id] <=> [b.updated_at, b.id] }
      first_file = processed_files.min { |a, b| [a.updated_at, a.id] <=> [b.updated_at, b.id] }
      ids = processed_files.map(&:id)

      if last_file.nil?
        return
      end

      cursor[:limit] = file_request.limit
      if cursor.include?(first_file)
        cursor[:ids] += ids
        cursor[:offset] = file_request.offset
      else
        cursor[:ids] = ids
        cursor[:limit] = nil
      end

      cursor.updated_at = last_file.updated_at
    end

    def execute(subcommand)
      if subcommand == 'reset'
        logger.warn("Resetting etl cursors!")
        @cursor_group.reset_all!
      elsif subcommand == 'run'
        while true
          return unless run_once
          logger.info("Continuing to process...")
        end

        logger.info("No more work detected, resting")
      else
        raise "Unknown subcommand #{subcommand}, did you mean #{completions_for('subcommand').join(' | ')}?"
      end
    end

    def completions_for(parameter)
      if parameter == 'subcommand'
        ['reset', 'run']
      else
        []
      end
    end

    def setup(config)
      super
      Polyphemus.instance.setup_db
    end
  end
end
