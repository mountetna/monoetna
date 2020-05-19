class FileController < Metis::Controller
  def remove
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::Forbidden, 'File is read-only' if file.read_only?

    raise Etna::Forbidden, 'Folder is read-only' if file.folder&.read_only?

    response = { files: [ file.to_hash ] }

    file.remove!

    return success_json(response)
  end

  def protect
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'File is already read-only' if file.read_only?

    file.protect!

    success_json(files: [ file.to_hash(@request) ])
  end

  def unprotect
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'File is not protected' unless file.read_only?

    file.unprotect!

    success_json(files: [ file.to_hash(@request) ])
  end

  def rename
    require_param(:new_file_path)
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::Forbidden, 'File is read-only' if file.read_only?

    raise Etna::Forbidden, 'File is read-only' if file.folder&.read_only?

    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:new_file_path])

    new_folder_path, new_file_name = Metis::File.path_parts(@params[:new_file_path])

    new_folder = require_folder(bucket, new_folder_path)

    raise Etna::Forbidden, 'Folder is read-only' if new_folder&.read_only?

    raise Etna::Forbidden, 'Cannot rename over existing file' if Metis::File.exists?(new_file_name, bucket, new_folder)

    raise Etna::Forbidden, 'Cannot rename over existing folder' if  Metis::Folder.exists?(new_file_name, bucket, new_folder)

    file.rename!(new_folder, new_file_name)

    success_json(files: [ file.to_hash(@request) ])
  end

  def copy
    require_param(:new_file_path)
    bucket = require_bucket

    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:new_file_path])

    raise Etna::Forbidden, "Cannot edit project #{@params[:project_name]}" unless @user.can_edit?(@params[:project_name])

    new_folder_path, new_file_name = Metis::File.path_parts(@params[:new_file_path])

    new_bucket = require_bucket(@params[:new_bucket_name])

    new_folder = require_folder(new_bucket, new_folder_path)

    raise Etna::Forbidden, 'Folder is read-only' if new_folder&.read_only?

    raise Etna::Forbidden, 'Cannot copy over existing file' if Metis::File.exists?(new_file_name, new_bucket, new_folder)

    raise Etna::Forbidden, 'Cannot copy over existing folder' if  Metis::Folder.exists?(new_file_name, new_bucket, new_folder)

    new_file = Metis::File.create(
      project_name: @params[:project_name],
      file_name: new_file_name,
      folder_id: new_folder&.id,
      bucket: new_bucket,
      author: Metis::File.author(@user),
      data_block: file.data_block
    )
    success_json(files: [ new_file.to_hash(@request) ])
  end

  def bulk_update
    # We want to support bulk copy / removal operations
    #   from Magma.
    # First, check that the user has access to all of
    #   the revision buckets and files, and that
    #   all revisions are valid.
    # If all revisions are valid, execute them.
    require_param(:revisions)

    revisions = @params[:revisions]

    raise Etna::BadRequest, 'At least one revision required' unless revisions.length > 0

    raise Etna::BadRequest, 'All revisions require "source" parameter' unless revisions.all? {|revision| revision.has_key? 'source'}
    raise Etna::BadRequest, 'All revisions require "dest" parameter' unless revisions.all? {|revision| revision.has_key? 'dest'}

    all_source_buckets = revisions.map {|rev| extract_bucket_from_path(rev[:source])}.uniq
    all_dest_buckets = revisions.map {|rev| extract_bucket_from_path(rev[:dest])}.uniq

    user_authorized_buckets = Metis::Bucket.where(
      project_name: @params[:project_name]
    ).all.select{|b| b.allowed?(@user, @request.env['etna.hmac'])}

    raise Etna::Forbidden, 'Cannot access the source buckets' unless all_source_buckets.all? {|bucket| user_authorized_buckets.include? bucket}
    raise Etna::Forbidden, 'Cannot access the destination buckets' unless all_dest_buckets.all? {|bucket| user_authorized_buckets.include? bucket}

    raise Etna::Forbidden, "Cannot edit project #{@params[:project_name]}" unless @user.can_edit?(@params[:project_name])

    # I'm not sure there is a way to check Files in bulk,
    #   without iterating through them all, like this?
    revisions.each do |revision|
      source_bucket = extract_bucket_from_path(revision[:source])

      file = Metis::File.from_path(
        source_bucket,
        extract_file_path_from_path(revision[:source]))

      raise Etna::Error.new('File not found', 404) unless file&.has_data?

      raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(
        extract_file_path_from_path(revision[:dest]))

      # Here we check for removal or linkage
      if revision[:dest].start_with? 'metis://'
        new_folder_path, new_file_name = Metis::File.path_parts(
          extract_file_path_from_path(revision[:dest]))

        new_bucket = require_bucket(extract_bucket_from_path(revision[:dest]))

        new_folder = require_folder(new_bucket, new_folder_path)

        raise Etna::Forbidden, '#{new_folder} Folder is read-only' if new_folder&.read_only?

        raise Etna::Forbidden, 'Cannot copy over existing file #{new_file_name}' if Metis::File.exists?(new_file_name, new_bucket, new_folder)

        raise Etna::Forbidden, 'Cannot copy over existing folder #{new_folder_path}' if  Metis::Folder.exists?(new_file_name, new_bucket, new_folder)
      else
        dest_bucket = extract_bucket_from_path(revision[:dest])
        dest_file = Metis::File.from_path(
          dest_bucket,
          extract_file_path_from_path(revision[:dest]))

        raise Etna::Error.new('Destiination file #{dest_file} not found', 404) unless dest_file&.has_data?

        raise Etna::Forbidden, 'Destination file #{dest_file} is read-only' if dest_file.read_only?

        raise Etna::Forbidden, 'Destination folder #{dest_file.folder} is read-only' if dest_file.folder&.read_only?
      end
    end

    # If we've gotten here, every revision looks good and we can execute them!
    new_files = []
    revisions.each do |revision|
      source_bucket = extract_bucket_from_path(revision[:source])

      file = Metis::File.from_path(
        source_bucket,
        extract_file_path_from_path(revision[:source]))

      # Here we check for removal or linkage
      if revision[:dest].start_with? 'metis://'
        new_folder_path, new_file_name = Metis::File.path_parts(
          extract_file_path_from_path(revision[:dest]))

        new_bucket = require_bucket(extract_bucket_from_path(revision[:dest]))

        new_folder = require_folder(new_bucket, new_folder_path)

        new_file = Metis::File.create(
          project_name: @params[:project_name],
          file_name: new_file_name,
          folder_id: new_folder&.id,
          bucket: new_bucket,
          author: Metis::File.author(@user),
          data_block: file.data_block
        )
        new_files.push(new_file.to_hash(@request))
      else
        dest_bucket = extract_bucket_from_path(revision[:dest])
        dest_file = Metis::File.from_path(
          dest_bucket,
          extract_file_path_from_path(revision[:dest]))

        new_files.push(dest_file.to_hash)

        dest_file.remove!
      end
    end
    success_json(files: new_files)
  end
end
