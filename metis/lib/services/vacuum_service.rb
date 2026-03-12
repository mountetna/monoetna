class Metis
  class VacuumService
    def initialize(project_name:, commit: false, user: nil)
      @project_name = project_name
      @commit = commit
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
          project_name: @project_name
        }
      }
    end

    private

    def find_orphaned_datablocks
      if @is_backfilled
        Metis::DataBlockLedger.find_orphaned_datablocks_backfilled_for_vacuum
      else
        Metis::DataBlockLedger.find_orphaned_datablocks_for_vacuum(@project_name)
      end
    end
  end
end
