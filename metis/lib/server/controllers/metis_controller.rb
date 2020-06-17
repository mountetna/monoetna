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

      raise Etna::BadRequest, "Invalid bucket: \"#{bucket_name || @params[:bucket_name]}\"" unless bucket

      raise Etna::Forbidden, "Cannot access bucket: \"#{bucket_name || @params[:bucket_name]}\"" unless bucket.allowed?(@user, @request.env['etna.hmac'])

      return bucket
    end

    def require_folder(bucket, folder_path)
      # if there is no folder_path, nothing is required
      return nil unless folder_path

      folder = Metis::Folder.from_path(bucket, folder_path).last

      raise Etna::BadRequest, "Invalid folder: \"#{folder_path}\"" unless folder

      return folder
    end
  end
end
