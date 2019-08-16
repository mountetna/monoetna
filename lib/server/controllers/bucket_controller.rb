class BucketController < Metis::Controller
  def list
    buckets = Metis::Bucket.where(
      project_name: @params[:project_name]
    ).all

    success_json(buckets: buckets.map(&:to_hash))
  end

  def create
    require_params(:owner, :access)

    raise Etna::BadRequest, 'Illegal bucket name' unless @params[:bucket_name] =~ /\A\w+\z/

    raise Etna::BadRequest, 'Invalid owner' unless Metis.instance.config(:hmac_keys).keys.include?(@params[:owner].to_sym)

    raise Etna::BadRequest, 'Invalid access' unless Metis::Bucket.valid_access?(@params[:access])

    bucket = Metis::Bucket.create(
      project_name: @params[:project_name],
      name: @params[:bucket_name],
      owner: @params[:owner],
      access: @params[:access]
    )
    success_json(bucket: bucket.to_hash)
  end

  def update
    bucket = require_bucket

    raise Etna::BadRequest, 'Invalid owner' unless Metis.instance.config(:hmac_keys).keys.include?(@params[:owner].to_sym)

    raise Etna::BadRequest, 'Invalid access' unless Metis::Bucket.valid_access?(@params[:access])


    bucket.owner = @params[:owner] if @params[:owner]
    bucket.access = @params[:access] if @params[:access]
    bucket.description = @params[:description] if @params[:description]
    bucket.save

    success_json(bucket: bucket.to_hash)
  end

  def remove
    bucket = require_bucket

    raise Etna::BadRequest, 'Cannot remove bucket' unless bucket.can_remove?

    bucket.remove!

    response = success_json(bucket: bucket.to_hash)
  end
end
