require_relative 'revision'

class Metis
  class FileRenameRevision < Revision
    def revise!
      raise 'Invalid revision, cannot revise!' unless valid?
      raise 'Cannot revise without a user' unless @user

      @source.file.update_bucket_and_rename!(
        @dest.folder,
        @dest.mpath.file_name,
        @dest.bucket,
        @user)
    end

    def validate_file(mpath_w_objs, file_check_type='source')
      errors_found = super

      # Do not let users clobber files by default
      file = Metis::File.from_folder(
        mpath_w_objs.bucket,
        mpath_w_objs.folder,
        mpath_w_objs.mpath.file_name)

      if file&.has_data? && file_check_type == 'dest'
        @errors.push("File \"#{mpath_w_objs.mpath.path}\" already exists")
        errors_found = true
      end

      if file&.read_only? && file_check_type == 'source'
        @errors.push("File \"#{mpath_w_objs.mpath.path}\" is read-only")
        errors_found = true
      end

      # Set this to use in revise!
      mpath_w_objs.file = file unless errors_found

      return errors_found
    end
  end
end