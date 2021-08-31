require_relative '../../query'

class BucketController < Metis::Controller
  def list
    buckets = Metis::Bucket.where(
      project_name: @params[:project_name]
    ).all.select{|b| b.allowed?(@user, @request.env['etna.hmac'])}

    success_json(buckets: buckets.map(&:to_hash))
  end

  def create
    require_params(:owner, :access)

    raise Etna::BadRequest, 'Illegal bucket name' unless Metis::Bucket.valid_bucket_name?(@params[:bucket_name])

    raise Etna::BadRequest, 'Invalid owner' unless Metis.instance.config(:hmac_keys).keys.include?(@params[:owner].to_sym)

    raise Etna::BadRequest, 'Invalid access' unless Metis::Bucket.valid_access?(@params[:access])

    raise Etna::BadRequest, 'Cannot create a reserved bucket' unless request_is_hmac_signed_by_owner?

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

  def find
    bucket = require_bucket
    require_params(:project_name, :params)
    params = @params[:params]

    limit = @params.has_key?(:limit) ? @params[:limit].to_i : nil
    offset = @params.has_key?(:offset) ? @params[:offset].to_i : nil

    raise Etna::BadRequest, "Invalid offset" if offset&.negative?
    raise Etna::BadRequest, "Invalid limit" if limit&.negative?

    query = Metis::Query.new(
      project_name: @params[:project_name],
      bucket: bucket,
      params: params,
      limit: limit,
      offset: offset,
    )

    results = query.execute

    files = results[:files]
    folders = results[:folders]

    file_hashes = @params[:hide_paths] ?
      file_hashes_without_paths(files: files) : 
      file_hashes_with_calculated_paths(
        files: files,
        bucket: bucket,
      )

    folder_hashes = @params[:hide_paths] ?
      folder_hashes_without_paths(folders: folders) :
      folder_hashes_with_calculated_paths(
        target_folders: folders,
        bucket: bucket
      )

    success_json(files: file_hashes, folders: folder_hashes)
  end

  private

  def request_is_hmac_signed_by_owner?
    return true unless @params[:owner].downcase == @params[:bucket_name].downcase

    return @hmac&.valid? && @hmac.id.to_s == @params[:owner]
  end
end
