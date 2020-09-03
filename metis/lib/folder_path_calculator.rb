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

    def get_folder_path(folder)
      # No parent folder, so just return this (root) folder's folder_name
      return folder.folder_name if !folder.folder_id

      # Use the cached path value for the parent folder if it exists
      return "#{@path_cache[folder.folder_id.to_s.to_sym]}/#{folder.folder_name}" if
        @path_cache.has_key?(folder.folder_id.to_s.to_sym)

      # Find the path for the parent folder, recursively
      parent_folder = @folders_by_id[folder.folder_id].first
      "#{get_folder_path(parent_folder)}/#{folder.folder_name}"
    end
  end
end