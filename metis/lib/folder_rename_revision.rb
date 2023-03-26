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

    def validate
      super
      validate_non_recursive_rename
    end

    private

    def validate_non_recursive_rename
      return unless dest_path_inside_source

      @errors << "Cannot move folder into itself: \"#{@source.mpath.path}\""
    end

    def dest_path_inside_source
      # Copying from root-to-root is okay
      return false if @source.folder.nil? && @dest.folder.nil?

      # Copying from some subfolder to root is okay
      return false if @dest.folder.nil?
      # Not okay where the source folder is a root folder and is the dest folder

      # Should return `true` if source is a folder, and it
      #   lies anywhere on the path of @dest.folder up to root,
      #   or source_folder == @dest.folder (which is parent of the final dest)
      @dest.folder == source_folder || @dest.folder.parent_folders.include?(source_folder)
    end

    def source_folder
      # source_folder is the @source.mpath.file_name as a Metis::Folder
      # If we get here past the other validations, there should
      #   be a Metis::Folder.
      query_params = {
        folder_name: @source.mpath.file_name,
        project_name: @source.mpath.project_name,
        bucket: @source.bucket
      }
      query_params[:folder] = @source.folder unless @source.folder.nil?

      Metis::Folder.where(*query_params).first
    end
  end
end
