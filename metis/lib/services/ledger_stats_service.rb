class Metis
  class LedgerStatsService
    def initialize(project_name: nil, backfilled: false, include_projects: [])
      @project_name = project_name
      @backfilled = backfilled
      @include_projects = Array(include_projects)
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
      blocked_by_project = build_blocked_by_project(blocked_datablocks)

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
      event_counts = Metis::DataBlockLedger.calculate_event_counts(@project_name, include_projects: @include_projects)

      # Planning view: candidates within the include_projects scope
      stats_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_for_stats(@project_name, include_projects: @include_projects)

      # Vacuum view: candidates safe to delete right now (global live-file check)
      ready_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_for_vacuum(@project_name)

      ready_ids = ready_datablocks.map(&:id).to_set
      blocked_datablocks = stats_datablocks.reject { |db| ready_ids.include?(db.id) }
      blocked_by_project = build_blocked_by_project(blocked_datablocks)

      project_breakdown = Metis::DataBlockLedger.calculate_project_breakdown(stats_datablocks, @project_name, @include_projects)

      vacuum_stats = {
        datablocks_ready: ready_datablocks.length,
        space_ready: ready_datablocks.sum(&:size),
        datablocks_blocked: blocked_datablocks.length,
        space_blocked: blocked_datablocks.sum(&:size),
        blocked_by_project: blocked_by_project,
        include_projects: @include_projects,
        project_breakdown: project_breakdown,
        details: Metis::DataBlockLedger.build_vacuum_details(stats_datablocks, @project_name)
      }

      {
        project_name: @project_name,
        event_counts: event_counts,
        vacuum: vacuum_stats
      }
    end

    def build_blocked_by_project(blocked_datablocks)
      return {} if blocked_datablocks.empty?

      blocked_by_project = {}
      Metis::File.where(data_block_id: blocked_datablocks.map(&:id)).each do |file|
        blocked_by_project[file.project_name] ||= 0
        blocked_by_project[file.project_name] += 1
      end
      blocked_by_project
    end
  end
end
