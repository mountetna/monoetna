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
    
    # Find orphaned datablocks for this project using the ledger
    orphaned_datablocks = Metis::DataBlockLedger.find_orphaned_datablocks(project_name)
    
    vacuumed = []
    errors = []
    
    orphaned_datablocks.each do |datablock|
      begin
        # Double-check it's not referenced (safety check)
        referenced_files = Metis::File.where(data_block_id: datablock.id).all
        if referenced_files.any?
          file_paths = referenced_files.map(&:file_path)
          errors << "DataBlock #{datablock.md5_hash} is still referenced by files: #{file_paths.join(', ')}"
          next
        end
        
        # Log vacuum event BEFORE deletion
        Metis::DataBlockLedger.log_vacuum(
          datablock: datablock,
          project_name: project_name,
          user: @user
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
      message: "vacuumed orphaned datablocks for project #{project_name}",
      payload: {
        vacuumed_count: vacuumed.length,
        space_freed: response[:summary][:space_freed]
      },
      consolidate: false
    )
    
    success_json(response)
  end
end
