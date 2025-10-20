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

  def delete_datablocks
    require_param(:md5_hashes)
    
    # Validate MD5 hash format
    raise Etna::BadRequest, "Invalid MD5 hash format!" unless @params[:md5_hashes].all? do |md5|
      md5 =~ /^[a-f0-9]{32}$/i
    end
    
    deleted_blocks = []
    errors = []
    
    @params[:md5_hashes].each do |md5_hash|
      begin
        # Find the datablock by MD5 hash
        datablock = Metis::DataBlock.where(md5_hash: md5_hash).first
        
        unless datablock
          errors << "DataBlock not found for MD5: #{md5_hash}"
          next
        end
        
        # Check if datablock is already removed
        if datablock.removed?
          errors << "DataBlock #{md5_hash} is already removed"
          next
        end
        
        # Check if datablock is still referenced by files
        referenced_files = Metis::File.where(data_block_id: datablock.id).all
        if referenced_files.any?
          file_paths = referenced_files.map(&:file_path)
          errors << "DataBlock #{md5_hash} is still referenced by files: #{file_paths.join(', ')}"
          next
        end
        
        # Delete the datablock
        datablock.remove!
        deleted_blocks << {
          md5_hash: datablock.md5_hash,
          size: datablock.actual_size,
          location: datablock.location
        }
        
      rescue => e
        errors << "Error deleting datablock #{md5_hash}: #{e.message}"
      end
    end
    
    response = {
      deleted: deleted_blocks,
      errors: errors,
      summary: {
        total_requested: @params[:md5_hashes].length,
        successfully_deleted: deleted_blocks.length,
        errors: errors.length
      }
    }
    
    success_json(response)
  end
end
