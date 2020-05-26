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

    def get_file_obj_from_path(path)
      # Assumes path is an instance of Metis::Path
      Metis::File.from_path(
        require_bucket(path.bucket_name),
        path.file_path)
    end

    def get_bucket_folder_file_from_path(path)
      # Assumes path is an instance of Metis::Path
      new_folder_path, new_file_name = Metis::File.path_parts(path.file_path)

      new_bucket = require_bucket(path.bucket_name)

      new_folder = require_folder(new_bucket, new_folder_path)

      return new_bucket, new_folder, new_file_name
    end
  end
end
