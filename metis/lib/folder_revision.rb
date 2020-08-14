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
      # Because of how the MetisPath regex works, the actual
      #   folder name gets returned as the mpath.file_name, so
      #   here we have to reconstruct the entire path by concatenating
      #   mpath.folder_path + mpath.file_name.
      mpath_w_objs.mpath.folder_path ?
        "#{mpath_w_objs.mpath.folder_path}/#{mpath_w_objs.mpath.file_name}" :
        mpath_w_objs.mpath.file_name
    end

    def validate_source_folder(mpath_w_objs)
      # NOTE: the "mpath_w_objs.folder" is the parent folder, not the
      #   actual folder, so we have to account for that in the validations and
      #   can't use the superclass #validate_folder()
      errors_found = false

      target_folder = Metis::Folder.from_path(
        mpath_w_objs.bucket,
        full_folder_path(mpath_w_objs)).last

      if !target_folder
        @errors.push("Folder not found: \"#{mpath_w_objs.mpath.path}\"")
        errors_found = true
      end

      if target_folder&.read_only?
        @errors.push("Folder \"#{mpath_w_objs.mpath.path}\" is read-only")
        errors_found = true
      end

      return errors_found
    end

    def validate_dest_is_not_file(mpath_w_objs)
      # Because we don't want users overwriting files with a folder.
      # Since we just need to check for file existence and don't
      #   care about the file.read_only? flag, we can't use
      #   the superclass #validate_file() method.
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
      return unless validate_source_folder(source_mpath_w_objs)
      return unless validate_folder(source_mpath_w_objs, 'source')
    end

    def validate_dest (dest_mpath_w_objs)
      return unless validate_mpath(dest_mpath_w_objs.mpath)
      return unless validate_bucket(dest_mpath_w_objs)
      return unless validate_folder(dest_mpath_w_objs, 'dest')
      return unless validate_dest_is_not_file(dest_mpath_w_objs)
    end
  end
end