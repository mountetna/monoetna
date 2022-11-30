require_relative 'folder_revision'

class Metis
  class FolderCopyRevision < FolderRevision
    def revise!
      raise 'Invalid revision, cannot revise!' unless valid?
      raise 'Cannot revise without a user' unless @user

      # Have to fetch the actual source folder. What gets
      #   set in the revision's source.folder attribute
      #   is the parent folder, which
      #   we don't need at this point.
      source_folder = Metis::Folder.from_path(
        @source.bucket, full_folder_path(@source)).last

      source_folder.copy(
        new_parent_folder: @dest.folder,
        new_bucket: @dest.bucket,
        user: @user)

      return source_folder
    end

    def validate
      super
      validate_non_recursive_copy
      validate_new_parent_folder
    end

    private

    def validate_new_parent_folder
      return if @source.bucket != @dest.bucket

      # Not okay if the source folder is a root folder and the dest is
      #   a root folder, too (i.e. dest.folder == nil), so would be
      #   copying into itself.
      return unless @source.folder != nil && @source.folder.folder_id == nil && @dest.folder == nil

      @errors << "Parent folder is the same as current: \"#{@source.mpath.path}\""
    end

    def validate_non_recursive_copy
      return unless dest_path_inside_source

      @errors << "Cannot copy folder into itself: \"#{@source.mpath.path}\""
    end

    def dest_path_inside_source
      # Copying from root-to-root is okay
      return false if @source.folder.nil? && @dest.folder.nil?

      # Otherwise, copying from some subfolder to root is okay
      return false if @dest.folder.nil?
      # Not okay where the source folder is a root folder and is the dest folder
      return @dest.folder.folder_name == @source.mpath.file_name && @dest.folder.folder.nil? if @source.folder.nil?

      return true if @source.folder == @dest.folder

      # Should return `true` if source is a folder, and it
      #   lies anywhere on the path of @dest.folder up to root
      @dest.folder.parent_folders.include?(@source.folder)
    end
  end
end