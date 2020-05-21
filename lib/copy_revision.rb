require_relative 'revision'

class Metis
    class CopyRevision < Revision
        def initialize(params)
            super(params)
            raise Etna::BadRequest, 'Copy revisions must have a non-nil "dest" parameter' unless params[:dest]
            raise Etna::BadRequest, 'Invalid dest path' unless valid_file_path?(params[:dest])
        end
    end
end