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

    raise Etna::Forbidden, 'Cannot copy over existing folder' if  Metis::Folder.exists?(new_file_name, new_bucket, new_folder)

    if Metis::File.exists?(new_file_name, new_bucket, new_folder)
      old_dest_file = Metis::File.from_path(new_bucket, @params[:new_file_path])
      old_dest_file.remove!
    end

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

  def bulk_copy
    # We want to support bulk copy operations
    #   from Magma.
    # First, check that the user has access to all of
    #   the revision buckets and files, and that
    #   all revisions are valid.
    # If all revisions are valid, execute them.
    require_param(:revisions)

    revisions = @params[:revisions]

    raise Etna::BadRequest, 'At least one revision required' unless revisions.length > 0

    raise Etna::BadRequest, 'All revisions require "source" parameter' unless revisions.all? {|revision| revision.has_key? :source}
    raise Etna::BadRequest, 'All revisions require "dest" parameter' unless revisions.all? {|revision| revision.has_key? :dest}

    all_source_buckets = revisions.select {|rev| rev[:source]}.
      map {|rev| extract_bucket_from_path(rev[:source])}.uniq
    all_dest_buckets = revisions.select {|rev| rev[:dest]}.
      map {|rev| extract_bucket_from_path(rev[:dest])}.uniq

    user_authorized_bucket_names = Metis::Bucket.where(
      project_name: @params[:project_name]
    ).all.select{|b| b.allowed?(@user, @request.env['etna.hmac'])}.
      map {|bucket| bucket.name}

    raise Etna::Forbidden, 'Cannot access the source buckets' unless all_source_buckets.all? {|bucket| user_authorized_bucket_names.include? bucket}
    raise Etna::Forbidden, 'Cannot access the destination buckets' unless all_dest_buckets.all? {|bucket| user_authorized_bucket_names.include? bucket}

    raise Etna::Forbidden, "Cannot edit project #{@params[:project_name]}" unless @user.can_edit?(@params[:project_name])

    revisions.each do |revision|
      file = get_file_obj_from_path(revision[:source])

      raise Etna::Error.new("File #{revision[:source]} not found", 404) unless file&.has_data?

      # Here we check for removal or linkage
      if revision[:dest]
        raise Etna::BadRequest, "Invalid path format for dest #{revision[:dest]}" unless revision[:dest].start_with? 'metis://'
        raise Etna::BadRequest, "Invalid path for dest #{revision[:dest]}" unless Metis::File.valid_file_path?(
          extract_file_path_from_path(revision[:dest]))

        new_bucket, new_folder, new_file_name = get_bucket_folder_file_from_path(revision[:dest])

        raise Etna::Forbidden, "#{new_folder.folder_name} folder is read-only" if new_folder&.read_only?

        # If the destination file exists, check if it is read-only, since
        #   we will remove it / modify it
        if Metis::File.exists?(new_file_name, new_bucket, new_folder)
          dest_file = get_file_obj_from_path(revision[:dest])
          raise Etna::Forbidden, "Destination file #{dest_file.file_name} is read-only" if dest_file.read_only?
        end

        raise Etna::Forbidden, "Cannot copy over existing folder #{revision[:dest]}" if  Metis::Folder.exists?(new_file_name, new_bucket, new_folder)
      end
    end

    # If we've gotten here, every revision looks good and we can execute them!
    new_files = []
    revisions.each do |revision|
      file = get_file_obj_from_path(revision[:source])

      # Here we check for removal or linkage
      if revision[:dest]
        new_bucket, new_folder, new_file_name = get_bucket_folder_file_from_path(revision[:dest])

        # If the destination file exists, remove it before creating
        #   the new link
        if Metis::File.exists?(new_file_name, new_bucket, new_folder)
          old_dest_file = get_file_obj_from_path(revision[:dest])
          old_dest_file.remove!
        end

        new_file = Metis::File.create(
          project_name: @params[:project_name],
          file_name: new_file_name,
          folder_id: new_folder&.id,
          bucket: new_bucket,
          author: Metis::File.author(@user),
          data_block: file.data_block
        )
        new_files.push(new_file.to_hash(@request))
      end
    end
    success_json(files: new_files)
  end
end
