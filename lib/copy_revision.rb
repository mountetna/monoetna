require_relative 'revision'

class Metis
    class CopyRevision < Revision
        def initialize(params)
            super(params)
            raise Etna::BadRequest, 'Copy revisions must have a non-nil "dest" parameter' unless params[:dest]
            raise Etna::BadRequest, "Invalid path for dest #{params[:dest]}" unless valid_file_path?(params[:dest])
        end

        def self.create_from_parts(params)
            Metis::CopyRevision.new({
              source: path_from_parts(
                params[:source][:project_name],
                params[:source][:bucket_name],
                params[:source][:file_path]),
              dest: path_from_parts(
                params[:dest][:project_name],
                params[:dest][:bucket_name],
                params[:dest][:file_path])
            })
          end
    end
end