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

  def tail
    bucket = require_bucket
    batch_start = @params[:batch_start]
    batch_end = @params[:batch_end]
    type = @params[:type]
    begin
      batch_start = DateTime.parse(batch_start)
      batch_end = DateTime.parse(batch_end)
    rescue
      raise Etna::BadRequest, 'Invalid batch end or start'
    end
    raise Etna::BadRequest, "Invalid type, must be one of 'folders' or 'files'" unless ['folders', 'files'].include?(type)

    try_stream('application/x-json-stream') do |stream|
      if type == 'files'
        head_query = <<QUERY
    SELECT 
      'file' as type, files.folder_id as parent_id, files.file_name as node_name, files.updated_at, files.created_at, data_blocks.md5_hash as file_hash, data_blocks.archive_id
    FROM files
    JOIN data_blocks ON data_blocks.id = files.data_block_id
    WHERE files.updated_at >= :start AND files.updated_at <= :end AND files.bucket_id = :bucket_id
QUERY
      else
        head_query = <<QUERY
    SELECT 
      'folder' as type, folder_id as parent_id, folder_name as node_name, updated_at, created_at, NULL as file_hash, NULL as archive_id
    FROM folders
    WHERE folders.updated_at >= :start AND folders.updated_at <= :end AND folders.bucket_id = :bucket_id
QUERY
      end

      results = Metis.instance.db[<<QUERY, {start: batch_start, end: batch_end, bucket_id: bucket.id}]
  WITH RECURSIVE tail_nodes AS (
#{head_query}
    UNION
    SELECT
      'parent' as type, folders.folder_id as parent_id, folders.folder_name as node_name, folders.updated_at, folders.created_at, NULL as file_hash, NULL as archive_id
    FROM folders
      JOIN tail_nodes ON folders.id = tail_nodes.parent_id
  )
  SELECT * FROM tail_nodes
QUERY
      results.stream.each do |row|
        JSON.dump(row, stream)
        stream.write("\n")
      end
    end
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

    file_hashes = hide_paths? ?
      file_hashes_without_paths(files: files) : 
      file_hashes_with_calculated_paths(
        files: files,
        bucket: bucket,
      )

    folder_hashes = hide_paths? ?
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

  def hide_paths?
    !!@params[:hide_paths]
  end
end
