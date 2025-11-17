class DataBlockController < Metis::Controller
  def exists
    raise Etna::BadRequest, "Improper md5!" unless @params[:md5s].all? do |md5|
      md5 =~ /^[a-f0-9]{32}$/i
    end
    if @params[:project_name]
      found = Metis::DataBlock.where(md5_hash: @params[:md5s]).join(
        Metis::File.where(project_name: @params[:project_name]),
        data_block_id: :id
      ).select_map(:md5_hash)
    else
      found = Metis::DataBlock.where(md5_hash: @params[:md5s]).select_map(:md5_hash)
    end
    success_json(found: found, missing: @params[:md5s] - found)
  end

  def vacuum_datablocks
    project_name = @params[:project_name]
    is_legacy = (project_name == 'legacy')
    
    # Find orphaned datablocks using the ledger
    orphaned_datablocks = if is_legacy
      Metis::DataBlockLedger.find_orphaned_datablocks_legacy
    else
      Metis::DataBlockLedger.find_orphaned_datablocks(project_name)
    end
    
    vacuumed = []
    errors = []
    
    orphaned_datablocks.each do |datablock|
      begin

        # Log vacuum event BEFORE deletion
        Metis::DataBlockLedger.log_vacuum(
          datablock,
          is_legacy ? nil : project_name,
          @user
        )
        
        # Remove the datablock
        datablock.remove!
        
        vacuumed << {
          md5_hash: datablock.md5_hash,
          size: datablock.size,
          location: datablock.location
        }
      rescue => e
        errors << "Failed to vacuum datablock #{datablock.md5_hash}: #{e.message}"
      end
    end
    
    response = {
      vacuumed: vacuumed,
      errors: errors,
      summary: {
        total_vacuumed: vacuumed.length,
        space_freed: vacuumed.sum { |d| d[:size] || 0 },
        errors_count: errors.length,
        project_name: project_name
      }
    }
    
    event_log(
      event: 'vacuum_datablocks',
      message: is_legacy ? "vacuumed legacy orphaned datablocks" : "vacuumed orphaned datablocks for project #{project_name}",
      payload: {
        vacuumed_count: vacuumed.length,
        space_freed: response[:summary][:space_freed]
      },
      consolidate: false
    )
    
    success_json(response)
  end
end
