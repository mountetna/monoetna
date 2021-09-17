require_relative "../../folder_path_calculator"

class Metis
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path("../views", __dir__)

    def success_json(obj)
      success(obj.to_json, "application/json")
    end

    def require_bucket(bucket_name = nil, project_name = nil)
      bucket = Metis::Bucket.where(
        project_name: project_name || @params[:project_name],
        name: bucket_name || @params[:bucket_name],
      ).first

      raise Etna::BadRequest, "Invalid bucket on project #{project_name || @params[:project_name]}: \"#{bucket_name || @params[:bucket_name]}\"" unless bucket

      raise Etna::Forbidden, "Cannot access bucket: \"#{bucket_name || @params[:bucket_name]}\"" unless bucket.allowed?(@user, @request.env["etna.hmac"])

      return bucket
    end

    def require_folder(bucket, folder_path)
      # if there is no folder_path, nothing is required
      return nil unless folder_path

      folder = Metis::Folder.from_path(bucket, folder_path).last

      raise Etna::BadRequest, "Invalid folder: \"#{folder_path}\"" unless folder

      return folder
    end

    def require_folder_by_id(bucket, folder_id)
      return nil unless folder_id

      begin
        folder_id = Integer(folder_id)
      rescue ArgumentError
        raise Etna::BadRequest, "folder_id must be numeric"
      end

      folder = Metis::Folder.where(id: folder_id).first

      raise Etna::BadRequest, "Invalid folder: \"#{folder_id}\"" unless folder

      return folder
    end

    private

    def file_hashes_without_paths(files:)
      return [] unless files

      files.map { |file|
        file.to_hash(with_path: false)
      }
    end

    def file_hashes_with_calculated_paths(bucket:, files:)
      # Calculate the folder_path, instead of
      #   doing it in the database because there is no bulk way to do that,
      #   and is recursive per file. Try to save queries / connections
      folder_path_calc = Metis::FolderPathCalculator.new(bucket: bucket)

      return [] unless files

      files.map { |file|
        if file.folder_id
          path = "#{folder_path_calc.get_folder_path(file.folder_id)}/#{file.file_name}"
        else
          path = file.file_name
        end

        file_hash = file.to_hash(file_path: path, request: @request)
        file_hash
      }
    end

    def folder_hashes_without_paths(folders:)
      return [] unless folders

      folders.map { |folder|
        folder.to_hash(false)
      }
    end

    def folder_hashes_with_calculated_paths(bucket:, limit: nil, offset: 0, all_folders: nil, target_folders: nil)
      # Calculate the folder_path, instead of
      #   doing it in the database.
      folder_path_calc = Metis::FolderPathCalculator.new(all_folders: all_folders, bucket: bucket)

      # Sorting folders by depth level makes some subsequent calculations simpler,
      #   especially when not paging. Shallow -> deep
      if all_folders
        sorted_folders = sort_folders_by_depth(
          limit: limit,
          offset: offset,
          all_folders: all_folders,
        )
        paged_folders = sorted_folders.slice(offset, limit)
      else
        paged_folders = target_folders
      end

      return [] unless paged_folders

      paged_folders.map { |fold|
        folder_hash = fold.to_hash(false)
        folder_hash[:folder_path] = folder_path_calc.get_folder_path(fold.id)
        folder_hash
      }
    end

    def sort_folders_by_depth(all_folders:, limit:, offset:)
      # Sorting folders by depth level makes some subsequent calculations simpler,
      #   especially when not paging. Shallow -> deep
      sorted_folders = []
      parent_folder_ids = [nil]
      folders_by_folder_id = all_folders.group_by { |fold| fold.folder_id }

      loop do
        child_folders = folders_by_folder_id.values_at(*parent_folder_ids).flatten.compact

        break if child_folders.length == 0

        parent_folder_ids = child_folders.map { |fold| fold.id }

        # Sort by folder_name within each depth level
        #   ... trying to make pagination consistent.
        sorted_folders += child_folders.sort { |f1, f2|
          f1[:folder_name] <=> f2[:folder_name]
        }

        break if sorted_folders.length >= limit + offset
      end
      sorted_folders
    end
  end
end
