require_relative 'revision'

class Metis
    class CopyRevision < Revision
        def self.create_from_parts(params)
            Metis::CopyRevision.new({
              source: Metis::Path.path_from_parts(
                params[:source][:project_name],
                params[:source][:bucket_name],
                params[:source][:file_path]),
              dest: Metis::Path.path_from_parts(
                params[:dest][:project_name],
                params[:dest][:bucket_name],
                params[:dest][:file_path])
            })
        end

        def valid? (validation_type=nil, user_authorized_bucket_names=nil)
            case validation_type
            when 'dest_path'
              return @dest.valid?
            when 'dest_bucket_access'
              return false unless user_authorized_bucket_names
              return user_authorized_bucket_names.include? @dest.bucket_name
            end
            super(validation_type, user_authorized_bucket_names)
        end
    end
end