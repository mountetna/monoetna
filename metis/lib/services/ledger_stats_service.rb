class Metis
  class LedgerStatsService
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

      vacuum_stats = {
        datablocks_ready: ready_datablocks.length,
        space_ready: ready_datablocks.sum(&:size),
        datablocks_blocked: blocked_datablocks.length,
        space_blocked: blocked_datablocks.sum(&:size),
        blocked_by_project: blocked_by_project,
        details: Metis::DataBlockLedger.build_vacuum_details(stats_datablocks)
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
      blocked_result = Metis::DataBlockLedger.find_blocked_datablocks(@project_name)
      blocked_datablocks = blocked_result[:datablocks]
      blocked_by = blocked_result[:blocked_by]

      vacuum_stats = {
        datablocks_ready: ready_datablocks.length,
        space_ready: ready_datablocks.sum(&:size),
        datablocks_blocked: blocked_datablocks.length,
        space_blocked: blocked_datablocks.sum(&:size),
        blocked_by_project: build_blocked_by_project_summary(blocked_datablocks, blocked_by),
        details: Metis::DataBlockLedger.build_vacuum_details(ready_datablocks, @project_name)
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
