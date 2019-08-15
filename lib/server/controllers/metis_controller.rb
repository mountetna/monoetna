class Metis
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def success_json(obj)
      success(obj.to_json, 'application/json')
    end

    def require_bucket(user_required=true)
      bucket = Metis::Bucket.where(
        project_name: @params[:project_name],
        name: @params[:bucket_name]
      ).first


      raise Etna::BadRequest, 'Invalid bucket!' unless bucket

      raise Etna::Forbidden, 'Cannot access bucket!' unless !user_required || bucket.allowed?(@user)

      return bucket
    end

    def require_folder(bucket, folder_path)
      # if there is no folder_path, nothing is required
      return nil unless folder_path

      folder = Metis::Folder.from_path(bucket, folder_path).last

      raise Etna::BadRequest, 'Invalid folder' unless folder

      return folder
    end
  end
end
