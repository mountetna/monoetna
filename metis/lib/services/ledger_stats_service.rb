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
      orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
      event_counts = Metis::DataBlockLedger.calculate_event_counts(nil)
      
      vacuum_stats = {
        datablocks_can_vacuum: orphaned_datablocks.length,
        space_can_clear: orphaned_datablocks.sum(&:size),
        details: Metis::DataBlockLedger.build_vacuum_details(orphaned_datablocks)
      }
      
      {
        project_name: 'backfilled',
        event_counts: event_counts,
        vacuum: vacuum_stats
      }
    end

    def calculate_tracked_stats
      event_counts = Metis::DataBlockLedger.calculate_event_counts(@project_name, include_projects: @include_projects)
      orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks(@project_name, include_projects: @include_projects)
      project_breakdown = Metis::DataBlockLedger.calculate_project_breakdown(orphaned_datablocks, @project_name, @include_projects)
      
      vacuum_stats = {
        datablocks_can_vacuum: orphaned_datablocks.length,
        space_can_clear: orphaned_datablocks.sum(&:size),
        include_projects: @include_projects,
        project_breakdown: project_breakdown,
        details: Metis::DataBlockLedger.build_vacuum_details(orphaned_datablocks, @project_name)
      }
      
      {
        project_name: @project_name,
        event_counts: event_counts,
        vacuum: vacuum_stats
      }
    end
  end
end

