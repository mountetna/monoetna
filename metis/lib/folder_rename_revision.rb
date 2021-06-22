require_relative 'folder_revision'

class Metis
  class FolderRenameRevision < FolderRevision
    def revise!
      raise 'Invalid revision, cannot revise!' unless valid?
      raise 'Cannot revise without a user' unless @user

      # Have to fetch the actual source folder. What gets
      #   set in the revision's source.folder attribute
      #   is the parent folder, which
      #   we don't need at this point.
      source_folder = Metis::Folder.from_path(
        @source.bucket, full_folder_path(@source)).last

      if @source.mpath.bucket_name == @dest.mpath.bucket_name
        source_folder.rename!(@dest.folder, @dest.mpath.file_name, @user)
      else
        source_folder.update_bucket_and_rename!(
          @dest.folder, @dest.mpath.file_name, @dest.bucket, @user)
      end

      return source_folder
    end
  end
end