require_relative 'revision'

class Metis
  class CopyRevision < Revision
    def initialize(params)
      super
      @restriction = params[:restriction]
    end

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

    def validate
      super
      if @restriction
        if @restriction == 'restricted' && !@source.data_block.restricted
          @errors.push("Unrestricted file \"#{@source.mpath.path}\" cannot be copied to a restricted destination")
          return
        end
        if @restriction == 'unrestricted' && @source.data_block.restricted
          @errors.push("Restricted file \"#{@source.mpath.path}\" cannot be copied to an unrestricted destination")
          return
        end
      end
    end
  end
end
