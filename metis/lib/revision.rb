require_relative 'path_with_objects'

class Metis
  class Revision
    attr_reader :source, :dest, :errors, :user
    def initialize(params)
      @source = Metis::PathWithObjects.new(params[:source])
      @dest = Metis::PathWithObjects.new(params[:dest])

      @paths = [@source, @dest]

      @user = params[:user]
      @errors = nil
      # TODO: Will have to handle DELETE differently,
      #       when params[:dest] will be nil
    end

    def set_bucket(path_with_objects, user_authorized_buckets)
      return unless path_with_objects.mpath.valid?
      path_with_objects.bucket = user_authorized_buckets.find do |b|
        path_with_objects.mpath.bucket_matches?(b)
      end
    end

    def set_folder(path_with_objects, bucket_folders)
      return unless path_with_objects.mpath.valid?
      path_with_objects.folder = bucket_folders.find do |f|
        path_with_objects.mpath.folder_matches?(f)
      end
    end

    def valid?
      !!@errors&.empty?
    end

    def validate
      @errors = []
      validate_source(@source)
      validate_dest(@dest)
    end

    def to_hash
      {
        dest: @dest.mpath.path,
        source: @source.mpath.path,
        errors: @errors
      }
    end

    def bucket_names
      # This is kind of weird, but we need the ability to grab
      #   all relevant bucket names, even before validation of
      #   the Revision.
      # So, if @source isn't valid, we return empty list.
      @paths.map do |path|
        path.mpath.bucket_name if path.mpath.valid?
      end.compact
    end

    def mpaths
      @paths.map do |path|
        path.mpath if path.mpath.valid?
      end.compact
    end

    private
    def validate_mpath(mpath)
      if mpath.valid?
        return true
      end
      @errors.push("Invalid path: \"#{mpath.path}\"")
      return false
    end

    def validate_bucket(mpath_w_objs)
      # Not having a bucket set means either it's a non-existent
      #   bucket or the user doesn't have access to it.
      if mpath_w_objs.bucket
        return true
      end
      @errors.push("Invalid bucket \"#{mpath_w_objs.mpath.bucket_name}\" in project \"#{mpath_w_objs.mpath.project_name}\". Check the bucket name and your permissions.")
      return false
    end

    def validate_file(mpath_w_objs, file_check_type='source')
      errors_found = false

      file = Metis::File.from_folder(
        mpath_w_objs.bucket,
        mpath_w_objs.folder,
        mpath_w_objs.mpath.file_name)

      if !file&.has_data? && file_check_type == 'source'
        @errors.push("File \"#{mpath_w_objs.mpath.path}\" not found")
        errors_found = true
      end

      if file&.restrict_user?(@user) && file_check_type == 'source'
        @errors.push("File \"#{mpath_w_objs.mpath.path}\" is restricted")
        errors_found = true
      end


      if file&.read_only? && file_check_type == 'dest'
        @errors.push("File \"#{mpath_w_objs.mpath.path}\" is read-only")
        errors_found = true
      end

      # Set this to use in revise!
      mpath_w_objs.file = file unless errors_found
      mpath_w_objs.data_block = file.data_block unless errors_found || !file

      return errors_found
    end

    def validate_folder(mpath_w_objs, folder_check_type='source')
      errors_found = false

      if (Metis::Folder.exists?(
          mpath_w_objs.mpath.file_name,
          mpath_w_objs.bucket,
          mpath_w_objs.folder) && folder_check_type == 'dest')

          @errors.push(
            "Cannot write over existing folder: \"#{mpath_w_objs.mpath.path}\""
          )
          errors_found = true
      end

      if mpath_w_objs.mpath.folder_path
        if !mpath_w_objs.folder
          @errors.push(
            "Invalid folder: \"#{mpath_w_objs.mpath.folder_path}\""
          )
          errors_found = true
        end
        if mpath_w_objs.folder&.read_only?
          @errors.push(
            "Folder \"#{mpath_w_objs.mpath.folder_path}\" is read-only"
          )
          errors_found = true
        end
      end

      return !errors_found
    end

    def validate_source (source_mpath_w_objs)
      return unless validate_mpath(source_mpath_w_objs.mpath)
      return unless validate_bucket(source_mpath_w_objs)
      return unless validate_file(source_mpath_w_objs, 'source')
    end

    def validate_dest (dest_mpath_w_objs)
      return unless validate_mpath(dest_mpath_w_objs.mpath)
      return unless validate_bucket(dest_mpath_w_objs)
      return unless validate_folder(dest_mpath_w_objs, 'dest')
      return unless validate_file(dest_mpath_w_objs, 'dest')
    end
  end
end
