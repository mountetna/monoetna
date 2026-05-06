class Metis
  class LedgerStatsService
    SIZE_DISTRIBUTION_BUCKETS = [
      { label: "0B-1MB", max: 1 * 1024 * 1024 },
      { label: "1MB-10MB", max: 10 * 1024 * 1024 },
      { label: "10MB-100MB", max: 100 * 1024 * 1024 },
      { label: "100MB-1GB", max: 1 * 1024 * 1024 * 1024 },
      { label: "1GB+", max: nil }
    ].freeze

    def initialize(project_name: nil, backfilled: false)
      @project_name = project_name
      @backfilled = backfilled
    end

    def calculate_stats
      if @backfilled
        calculate_backfilled_stats
      elsif @project_name
        calculate_tracked_stats
      else
        raise Etna::BadRequest, "Must provide either 'project_name' (for tracked mode) or 'backfilled' parameter"
      end
    end

    private

    def calculate_backfilled_stats
      stats_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled_for_stats
      ready_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled_for_vacuum
      event_counts = Metis::DataBlockLedger.calculate_event_counts(nil)

      ready_ids = ready_datablocks.map(&:id).to_set
      blocked_datablocks = stats_datablocks.reject { |db| ready_ids.include?(db.id) }
      blocked_by_project = backfilled_blocked_by_project(blocked_datablocks)
      full_metadata = build_backfilled_metadata(stats_datablocks)

      vacuum_stats = {
        datablocks_ready: ready_datablocks.length,
        space_ready: ready_datablocks.sum(&:size),
        datablocks_blocked: blocked_datablocks.length,
        space_blocked: blocked_datablocks.sum(&:size),
        blocked_by_project: blocked_by_project,
        full_metadata: full_metadata,
        date_distribution: build_date_distribution(full_metadata),
        size_distribution: build_size_distribution(full_metadata),
        extension_distribution: build_extension_distribution(full_metadata)
      }

      {
        project_name: 'backfilled',
        event_counts: event_counts,
        vacuum: vacuum_stats
      }
    end

    def calculate_tracked_stats
      event_counts = Metis::DataBlockLedger.calculate_event_counts(@project_name)

      ready_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_for_vacuum(@project_name)
      blocked_result = Metis::DataBlockLedger.find_blocked_orphaned_datablocks(@project_name)
      blocked_datablocks = blocked_result[:datablocks]
      blocked_by = blocked_result[:blocked_by]

      vacuum_stats = {
        datablocks_ready: ready_datablocks.length,
        space_ready: ready_datablocks.sum(&:size),
        datablocks_blocked: blocked_datablocks.length,
        space_blocked: blocked_datablocks.sum(&:size),
        blocked_by_project: build_blocked_by_project_summary(blocked_datablocks, blocked_by),
        full_metadata: Metis::DataBlockLedger.build_vacuum_details(ready_datablocks, @project_name)
      }

      {
        project_name: @project_name,
        event_counts: event_counts,
        vacuum: vacuum_stats
      }
    end

    def backfilled_blocked_by_project(blocked_datablocks)
      return {} if blocked_datablocks.empty?

      summary = {}
      Metis::File.where(data_block_id: blocked_datablocks.map(&:id)).each do |file|
        summary[file.project_name] ||= 0
        summary[file.project_name] += 1
      end
      summary
    end

    def build_backfilled_metadata(datablocks)
      datablock_ids = datablocks.map(&:id)
      return [] if datablock_ids.empty?

      events_by_id = Metis::DataBlockLedger
        .where(data_block_id: datablock_ids)
        .all
        .group_by(&:data_block_id)

      datablocks.map do |datablock|
        files = build_files_from_events(events_by_id[datablock.id] || [])

        {
          data_block_id: datablock.id,
          md5_hash: datablock.md5_hash,
          size: datablock.size,
          created_at: datablock.created_at&.iso8601,
          description: datablock.description,
          files: files,
          file_types: extract_extensions(files: files, description: datablock.description)
        }
      end
    end

    def build_files_from_events(events)
      events.map do |event|
        next if event.file_path.nil? && event.bucket_name.nil?

        {
          file_path: event.file_path,
          bucket_name: event.bucket_name
        }
      end.compact.uniq
    end

    def build_date_distribution(full_metadata)
      full_metadata.each_with_object(Hash.new(0)) do |record, distribution|
        key = record[:created_at]&.split("T")&.first || "unknown"
        distribution[key] += 1
      end
    end

    def build_size_distribution(full_metadata)
      distribution = SIZE_DISTRIBUTION_BUCKETS.each_with_object({}) { |bucket, h| h[bucket[:label]] = 0 }

      full_metadata.each do |record|
        size = record[:size].to_i
        bucket = SIZE_DISTRIBUTION_BUCKETS.find { |b| b[:max].nil? || size < b[:max] }
        distribution[bucket[:label]] += 1
      end

      distribution
    end

    def build_extension_distribution(full_metadata)
      full_metadata.each_with_object(Hash.new(0)) do |record, distribution|
        extensions = record[:file_types] || []
        if extensions.empty?
          distribution["unknown"] += 1
        else
          extensions.each { |ext| distribution[ext] += 1 }
        end
      end
    end

    def extract_extensions(files:, description:)
      extensions = []

      files.each do |file|
        ext = normalize_extension(file[:file_path])
        extensions << ext if ext
      end

      if extensions.empty?
        described_file = description&.match(/Originally for (.+)\z/)&.captures&.first
        ext = normalize_extension(described_file)
        extensions << ext if ext
      end

      extensions.uniq
    end

    def normalize_extension(file_name)
      return nil if file_name.nil? || file_name.empty?

      ext = ::File.extname(file_name).downcase.delete_prefix(".")
      ext.empty? ? "no_extension" : ext
    end

    # Summarizes blocked_by as { project_name => count_of_datablocks_it_blocks }
    def build_blocked_by_project_summary(blocked_datablocks, blocked_by)
      return {} if blocked_datablocks.empty?

      summary = {}
      blocked_datablocks.each do |db|
        (blocked_by[db.id] || []).each do |project|
          summary[project] ||= 0
          summary[project] += 1
        end
      end
      summary
    end
  end
end
