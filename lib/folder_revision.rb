require_relative 'revision'

class Metis
  class FolderRevision < Revision
    def validate
      @errors = []
      validate_source(@source)
      validate_dest(@dest)
    end

    def bucket_names
      # This is kind of weird, but we need the ability to grab
      #   all relevant bucket names, even before validation of
      #   the CopyRevision.
      # So, if @dest doesn't exist, we don't return it
      source_bucket_names = super
      if @dest.mpath.valid?
        return source_bucket_names.push(@dest.mpath.bucket_name)
      end
      return source_bucket_names
    end

    def mpaths
      source_mpaths = super
      source_mpaths.push(@dest.mpath) if @dest.mpath.valid?
      return source_mpaths
    end

    def to_hash
      {
        dest: @dest.mpath.path,
        source: @source.mpath.path,
        errors: @errors
      }
    end

    private

    def full_folder_path(mpath_w_objs)
      mpath_w_objs.mpath.folder_path ?
        "#{mpath_w_objs.mpath.folder_path}/#{mpath_w_objs.mpath.file_name}" :
        mpath_w_objs.mpath.file_name
    end

    def validate_target_folder(mpath_w_objs, folder_check_type='source')
      # NOTE: the "mpath_w_objs.folder" is the parent folder, not the
      #   actual folder, so we have to account for that in the validations.
      errors_found = false

      target_folder = Metis::Folder.from_path(
        mpath_w_objs.bucket,
        full_folder_path(mpath_w_objs)).last

      if !target_folder
        @errors.push("Folder not found: \"#{mpath_w_objs.mpath.path}\"")
        errors_found = true
      end

      if target_folder&.read_only? && folder_check_type == 'source'
        @errors.push("Folder \"#{mpath_w_objs.mpath.path}\" is read-only")
        errors_found = true
      end

      return errors_found
    end

    def validate_dest_file(mpath_w_objs)
      # Because we don't want users overwriting files with a folder
      errors_found = false

      file = Metis::File.from_folder(
        mpath_w_objs.bucket,
        mpath_w_objs.folder,
        mpath_w_objs.mpath.file_name)

      if file&.has_data?
        @errors.push("Cannot overwrite existing file: \"#{mpath_w_objs.mpath.path}\"")
        errors_found = true
      end

      return errors_found
    end

    def validate_source (source_mpath_w_objs)
      return unless validate_mpath(source_mpath_w_objs.mpath)
      return unless validate_bucket(source_mpath_w_objs)
      return unless validate_target_folder(source_mpath_w_objs, 'source')
      return unless validate_folder(source_mpath_w_objs, 'source')
    end

    def validate_dest (dest_mpath_w_objs)
      return unless validate_mpath(dest_mpath_w_objs.mpath)
      return unless validate_bucket(dest_mpath_w_objs)
      return unless validate_folder(dest_mpath_w_objs, 'dest')
      return unless validate_dest_file(dest_mpath_w_objs)
    end
  end
end