class Metis
  class Controller < Etna::Controller
    VIEW_PATH=File.expand_path('../views', __dir__)

    def success_json(obj)
      success(obj.to_json, 'application/json')
    end

    def require_bucket
      bucket = Metis::Bucket.where(
        project_name: @params[:project_name], 
        name: @params[:bucket_name]
      ).first

      raise Etna::BadRequest, 'Invalid bucket!' unless bucket && bucket.allowed?(@user)

      return bucket
    end
  end
end
