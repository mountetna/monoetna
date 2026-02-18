class Metis
  class VacuumService
    def initialize(project_name:, commit: false, include_projects: [], user: nil)
      @project_name = project_name
      @commit = commit
      @include_projects = Array(include_projects)
      @user = user
      @is_backfilled = (project_name == 'backfilled')
    end

    def vacuum_datablocks
      orphaned_datablocks = find_orphaned_datablocks
      vacuumed = []
      errors = []
      
      orphaned_datablocks.each do |datablock|
        begin
          if @commit
            datablock.remove!
            Metis::DataBlockLedger.log_vacuum(
              datablock,
              @is_backfilled ? nil : @project_name,
              @user
            )
          end
          
          vacuumed << {
            md5_hash: datablock.md5_hash,
            size: datablock.size,
            location: datablock.location
          }
        rescue => e
          errors << "Failed to vacuum datablock #{datablock.md5_hash}: #{e.message}"
        end
      end
      
      {
        dry_run: !@commit,
        vacuumed: vacuumed,
        errors: errors,
        summary: {
          total_vacuumed: vacuumed.length,
          space_freed: vacuumed.sum { |d| d[:size] || 0 },
          errors_count: errors.length,
          project_name: @project_name,
          include_projects: @include_projects
        }
      }
    end

    private

    def find_orphaned_datablocks
      if @is_backfilled
        Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
      else
        Metis::DataBlockLedger.find_orphaned_datablocks(@project_name, include_projects: @include_projects)
      end
    end
  end
end

