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

  def batch_query
    batch_start = @params[:batch_start]
    batch_end = @params[:batch_end]

    return nil if batch_start.nil? or batch_end.nil?

    begin
      batch_start = DateTime.parse(batch_start)
      batch_end = DateTime.parse(batch_end)
    rescue
      raise Etna::BadRequest, 'Invalid batch end or start'
    end

    ["head.updated_at >= :start AND head.updated_at <= :end", { start: batch_start, end: batch_end }]
  end

  def folder_query
    folder_id = @params[:folder_id]

    if folder_id.nil?
      ["head.folder_id IS NULL", {}]
    else
      ["head.folder_id = :folder_id", {folder_id: folder_id.to_i}]
    end
  end

  def tail_query
    batch_query || folder_query
  end

  def tail
    bucket = require_bucket
    type = @params[:type]
    raise Etna::BadRequest, "Invalid type, must be one of 'folders' or 'files'" unless ['folders', 'files'].include?(type)

    head_cond, head_params = tail_query

    try_stream('application/x-json-stream') do |stream|
      if type == 'files'
        head_query = <<QUERY
    SELECT 
      'file' as type, head.id as id, head.folder_id as parent_id, head.file_name as node_name, head.updated_at, head.created_at, data_blocks.md5_hash as file_hash, data_blocks.archive_id
    FROM files as head
    JOIN data_blocks ON data_blocks.id = head.data_block_id
    WHERE #{head_cond} AND head.bucket_id = :bucket_id
QUERY
      else
        head_query = <<QUERY
    SELECT 
      'folder' as type, head.id as id, folder_id as parent_id, folder_name as node_name, updated_at, created_at, NULL as file_hash, NULL as archive_id
    FROM folders as head
    WHERE #{head_cond} AND head.bucket_id = :bucket_id
QUERY
      end

      results = Metis.instance.db[<<QUERY, {bucket_id: bucket.id}.update(head_params)]
  WITH RECURSIVE tail_nodes AS (
#{head_query}
    UNION
    SELECT
      'parent' as type, folders.id as id, folders.folder_id as parent_id, folders.folder_name as node_name, folders.updated_at, folders.created_at, NULL as file_hash, NULL as archive_id
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
