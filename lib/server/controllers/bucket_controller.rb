class BucketController < Metis::Controller
  def list
    buckets = Metis::Bucket.where(
      project_name: @params[:project_name]
    ).all.select{|b| b.allowed?(@user)}

    success_json(buckets: buckets.map(&:to_hash))
  end

  def create
    require_params(:owner, :access)

    raise Etna::BadRequest, 'Illegal bucket name' unless Metis::Bucket.valid_bucket_name?(@params[:bucket_name])

    raise Etna::BadRequest, 'Invalid owner' unless Metis.instance.config(:hmac_keys).keys.include?(@params[:owner].to_sym)

    raise Etna::BadRequest, 'Invalid access' unless Metis::Bucket.valid_access?(@params[:access])

    existing_bucket = Metis::Bucket.where(
      project_name: @params[:project_name],
      name: @params[:bucket_name]
    ).first

    raise Etna::BadRequest, 'Duplicate bucket name' if existing_bucket

    bucket = Metis::Bucket.create(
      project_name: @params[:project_name],
      name: @params[:bucket_name],
      owner: @params[:owner],
      description: @params[:description],
      access: @params[:access]
    )
    bucket.create_actual_bucket!
    success_json(bucket: bucket.to_hash)
  end

  def update
    bucket = require_bucket

    raise Etna::BadRequest, 'Invalid access' if @params[:access] && !Metis::Bucket.valid_access?(@params[:access])

    raise Etna::BadRequest, 'Illegal bucket name' if @params[:new_bucket_name] && !Metis::Bucket.valid_bucket_name?(@params[:new_bucket_name])

    bucket.access = @params[:access] if @params[:access]
    bucket.description = @params[:description] if @params[:description]
    bucket.rename!(@params[:new_bucket_name]) if @params[:new_bucket_name]

    bucket.save

    success_json(bucket: bucket.to_hash)
  end

  def remove
    bucket = require_bucket

    raise Etna::BadRequest, 'Cannot remove bucket' unless bucket.can_remove?

    response = { bucket: bucket.to_hash }

    bucket.remove!

    return success_json(response)
  end
end
