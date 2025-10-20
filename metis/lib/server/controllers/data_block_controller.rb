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

  def get_datablocks
    bucket = require_bucket
    folder_path = @params[:folder_path]
    
    # Find the specific folder by path
    folder = Metis::Folder.from_path(bucket, folder_path).last
    return success_json(datablocks: [], message: "Folder not found") unless folder
    
    # Get all child folders recursively
    all_folders = [folder] + folder.child_folders
    folder_ids = all_folders.map(&:id)
    
    # Get all data blocks for files in all these folders
    results = Metis::DataBlock.left_join(:files, data_block_id: :id)
                             .where(Sequel[:files][:folder_id] => folder_ids)
                             .select(Sequel[:data_blocks][:md5_hash], 
                                    Sequel[:data_blocks][:size], 
                                    Sequel[:data_blocks][:description], 
                                    Sequel[:data_blocks][:created_at], 
                                    Sequel[:data_blocks][:removed],
                                    Sequel[:files][:file_name],
                                    Sequel[:files][:folder_id])
                             .all
    
    # Create the output data with full paths
    datablocks = results.map do |result|
      # Find the folder for this file
      file_folder = all_folders.find { |f| f.id == result[:folder_id] }
      
      # Build the full path: bucket_name/folder_path/file_name
      full_path = if file_folder
        folder_path_parts = file_folder.folder_path
        File.join(bucket.name, *folder_path_parts, result[:file_name])
      else
        File.join(bucket.name, result[:file_name])
      end
      
      {
        full_path: full_path,
        file_name: result[:file_name],
        md5_hash: result[:md5_hash],
        size: result[:size],
        description: result[:description],
        created_at: result[:created_at],
        removed: result[:removed]
      }
    end
    
    response = {
      datablocks: datablocks,
      summary: {
        total_datablocks: datablocks.length,
        folders_searched: all_folders.length,
        bucket_name: bucket.name,
        folder_path: folder_path
      }
    }
    
    success_json(response)
  end
end
