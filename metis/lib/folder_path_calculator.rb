# This helper
#   class calculates the paths in memory
#   instead of relying on database queries.

class Metis
  class FolderPathCalculator

    def initialize(all_folders: nil, bucket:)
      all_folders = all_folders ? all_folders : Metis::Folder.where(
        bucket: bucket
      ).all

      @folders_by_id = all_folders.group_by { |fold| fold.id }

      # To prevent too much redundant calculation,
      #   cache the path for discovered folders.
      @path_cache = {}
    end

    def get_folder_path(folder_id)
      folder = @folders_by_id[folder_id].first
      path = [folder.folder_name]
      folder_to_check = folder
      loop do
        # No parent folder, so just return this (root) folder's folder_name
        break if !folder_to_check.folder_id

        # Use the cached path value for the parent folder if it exists
        if @path_cache.has_key?(folder_to_check.folder_id.to_s.to_sym)
          path.unshift(*@path_cache[folder_to_check.folder_id.to_s.to_sym])
          break
        end

        # Add the parent folder to the path
        parent_folder = @folders_by_id[folder_to_check.folder_id].first
        path.unshift(parent_folder.folder_name)
        folder_to_check = parent_folder
      end
      # Once we've determined the path to a folder, we'll add it to our cache
      @path_cache[folder_id.to_s.to_sym] = path.dup
      path.join('/')
    end
  end
end