require_relative 'revision'

class Metis
  class CopyRevision < Revision
    def revise!
      raise 'Invalid revision, cannot revise!' unless valid?
      raise 'Cannot revise without a user' unless @user

      Metis::File.copy({
        project_name: @source.mpath.project_name,
        source_file: @source.file,
        dest_file_path: @dest.mpath.file_path,
        dest_bucket_name: @dest.mpath.bucket_name,
        user: @user
      })
    end
  end
end