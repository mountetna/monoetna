require_relative 'path_with_objects'

class Metis
  class Revision
    attr_reader :source, :dest, :errors, :user
    def initialize(params)
      @source = Metis::PathWithObjects.new(params[:source])
      @dest = Metis::PathWithObjects.new(params[:dest])

      @user = params[:user]
      @errors = nil
      # TODO: Will have to handle DELETE differently,
      #       when params[:dest] will be nil
    end

    def set_bucket(path_with_objects, user_authorized_buckets)
      path_with_objects.bucket = user_authorized_buckets.find {
        |b| b.name == path_with_objects.mpath.bucket_name
      }
    end

    def set_folder(path_with_objects, bucket_folders)
      path_with_objects.folder = bucket_folders.find {
        |f| f.folder_path == path_with_objects.mpath.folder_path
      }
    end

    def valid?
      @errors&.empty?
    end

    def validate (user_authorized_bucket_names)
      raise 'not implemented'
    end

    def bucket_names
      # This is kind of weird, but we need the ability to grab
      #   all relevant bucket names, even before validation of
      #   the Revision.
      # So, if @source isn't valid, we return empty list.
      if @source.mpath.valid?
        return [@source.mpath.bucket_name]
      end

      return []
    end

    def paths
      [@source.mpath]
    end

    private
    def validate_mpath(mpath)
      if mpath.valid?
        return true
      end
      @errors.push("Invalid path: #{mpath.path}")
      return false
    end

    def validate_bucket(mpath_w_objs)
      # Not having a bucket set means either it's a non-existent
      #   bucket or the user doesn't have access to it.
      if mpath_w_objs.bucket
        return true
      end
      @errors.push("Invalid bucket: #{mpath_w_objs.mpath.bucket_name}")
      return false
    end

    def validate_file(mpath_w_objs, check_read_only=false)
      errors_found = false

      file = Metis::File.from_folder(
        mpath_w_objs.bucket,
        mpath_w_objs.folder,
        mpath_w_objs.mpath.file_name)

      if !file&.has_data?
        @errors.push("File #{mpath_w_objs.mpath.path} not found")
        errors_found = true
      end

      if check_read_only && file&.read_only?
        @errors.push("File #{mpath_w_objs.mpath.path} is read-only")
        errors_found = true
      end

      # Set this to use in revise!
      mpath_w_objs.file = file unless errors_found

      return errors_found
    end

    def validate_folder(mpath_w_objs)
      errors_found = false

      if Metis::Folder.exists?(
          mpath_w_objs.mpath.file_name,
          mpath_w_objs.bucket,
          mpath_w_objs.folder)

          @errors.push(
            "Cannot copy over existing folder #{dest_mpath.path}"
          )
          errors_found = true
      end

      if mpath_w_objs.mpath.folder_path
        if !mpath_w_objs.folder
          @errors.push(
            "Invalid dest folder: #{mpath_w_objs.mpath.folder_path}"
          )
          errors_found = true
        end
        if mpath_w_objs.folder&.read_only?
          @errors.push(
            "Dest folder #{mpath_w_objs.mpath.folder_name} is read-only"
          )
          errors_found = true
        end
      end

      return !errors_found
    end

    def validate_source (source_mpath_w_objs)
      # First just check syntax of the path
      # Next check the bucket -- does it exist? Have access?
      #   When bulk setting buckets, if user not authorized to access the bucket or
      #   the bucket doesn't exist,
      #   just set it as nil. Then in validation check, any nil bucket
      #   returns "Invalid bucket"
      # Get all folders from the consolidated paths and buckets
      # Set these in the revisions / paths.
      # Then check the parent folder (dest only):
      #   dest: If it exists, if it is read-only
      # Then check the file_path.
      #   For source: does it exist and have data?
      #   For dest: is it read-only, is it really a folder?, etc.
      return unless validate_mpath(source_mpath_w_objs.mpath)
      return unless validate_bucket(source_mpath_w_objs)
      return unless validate_file(source_mpath_w_objs)
    end

    def validate_dest (dest_mpath_w_objs)
      return unless validate_mpath(dest_mpath_w_objs.mpath)
      return unless validate_bucket(dest_mpath_w_objs)
      return unless validate_folder(dest_mpath_w_objs)
      return unless validate_file(dest_mpath_w_objs, true)
    end
  end
end
