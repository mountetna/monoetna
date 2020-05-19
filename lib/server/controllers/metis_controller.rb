class Metis
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def success_json(obj)
      success(obj.to_json, 'application/json')
    end

    def require_bucket(bucket_name=nil)
      bucket = Metis::Bucket.where(
        project_name: @params[:project_name],
        name: bucket_name || @params[:bucket_name]
      ).first

      raise Etna::BadRequest, 'Invalid bucket' unless bucket

      raise Etna::Forbidden, 'Cannot access bucket' unless bucket.allowed?(@user, @request.env['etna.hmac'])

      return bucket
    end

    def require_folder(bucket, folder_path)
      # if there is no folder_path, nothing is required
      return nil unless folder_path

      folder = Metis::Folder.from_path(bucket, folder_path).last

      raise Etna::BadRequest, 'Invalid folder' unless folder

      return folder
    end

    def extract_bucket_from_path(path)
      # Assumes Metis file paths are in the form
      #    metis://<project>/<bucket>/<folder path>/<file name>
      # Splitting the above produces
      #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]
      # Should this be in a central gem, like etna, so
      #   we can share it across applications?
      metis_file_location_parts = revision[attribute_name].split('/')
      return revision[3]
    end

    def extract_file_path_from_path(path)
      # Assumes Metis file paths are in the form
      #    metis://<project>/<bucket>/<folder path>/<file name>
      # Splitting the above produces
      #   ["metis", "", "<project>", "<bucket>", "<folder path>" ... "file name"]
      # Should this be in a central gem, like etna, so
      #   we can share it across applications?
      metis_file_location_parts = revision[attribute_name].split('/')
      return metis_file_location_parts[4..-1].join('/'))
    end
  end
end
