require_relative 'revision'

class Metis
  class CopyRevision < Revision
    def revise!
      raise 'Invalid revision, cannot revise!' unless valid?
      raise 'Cannot revise without a user' unless @user

      Metis::File.link_to_block({
        project_name: @source.mpath.project_name,
        source_data_block: @source.data_block,
        dest_project_name: @dest.mpath.project_name,
        dest_file_path: @dest.mpath.file_path,
        dest_bucket_name: @dest.mpath.bucket_name,
        user: @user
      })
    end
  end
end