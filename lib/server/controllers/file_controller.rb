require_relative '../../copy_revision'

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

    # Keep this here to verify the bucket is valid,
    # since source bucket name is part of the path
    bucket = require_bucket

    raise Etna::Forbidden, "Cannot edit project #{@params[:project_name]}" unless @user.can_edit?(@params[:project_name])

    revision = Metis::CopyRevision.create_from_parts({
      source: {
        project_name: @params[:project_name],
        bucket_name: @params[:bucket_name],
        file_path: @params[:file_path]
      },
      dest: {
        project_name: @params[:project_name],
        bucket_name: @params[:new_bucket_name] ? @params[:new_bucket_name] : @params[:bucket_name],
        file_path: @params[:new_file_path]
      }
    })

    validate_copy_revision(revision)

    new_file = Metis::File.copy({
      project_name: @params[:project_name],
      source_file: revision.source_file,
      dest_file_path: revision.dest_file_path,
      dest_bucket_name: revision.dest_bucket_name,
      user: @user
    })

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

    raise Etna::Forbidden, "Cannot edit project #{@params[:project_name]}" unless @user.can_edit?(@params[:project_name])

    revisions = []

    JSON.parse(@params[:revisions]).
      each { |rev|
        revisions.push(Metis::CopyRevision.new(rev.transform_keys(&:to_sym)))
      }

    raise Etna::BadRequest, 'At least one revision required' unless revisions.length > 0

    # In bulk copy mode, we validate the buckets in bulk
    #   since they aren't part of the path (are part of revisions),
    #   and trying to minimize database hits by rejecting early, if
    #   possible.
    all_source_bucket_names = revisions.map {|rev| rev.source_bucket_name}.uniq
    all_dest_bucket_names = revisions.map {|rev| rev.dest_bucket_name}.uniq

    hmac = @request.env['etna.hmac']

    raise Etna::Forbidden 'Invalid signature' unless hmac.valid?

    user_authorized_bucket_names = Metis::Bucket.where(
      project_name: @params[:project_name],
      owner: ['metis', hmac.id.to_s],
      name: all_source_bucket_names + all_dest_bucket_names
    ).all.select{|b| b.allowed?(@user, hmac)}.
      map {|bucket| bucket.name}

    revisions.map {|rev| rev.validate_access_to_buckets(user_authorized_bucket_names)}

    revisions.each do |revision|
      validate_copy_revision(revision)
    end

    # If we've gotten here, every revision looks good and we can execute them!
    new_files = []
    revisions.each do |revision|

      new_file = Metis::File.copy({
        project_name: @params[:project_name],
        source_file: revision.source_file,
        dest_file_path: revision.dest_file_path,
        dest_bucket_name: revision.dest_bucket_name,
        user: @user
      })

      new_files.push(new_file.to_hash(@request))
    end
    success_json(files: new_files)
  end

  private

  def validate_copy_revision(revision)
    file = revision.source_file

    raise Etna::Error.new("File #{revision.source} not found", 404) unless file&.has_data?

    new_bucket, new_folder, new_file_name = get_bucket_folder_file_from_path(revision.dest)

    raise Etna::Forbidden, "#{new_folder.folder_name} folder is read-only" if new_folder&.read_only?

    # If the destination file exists, check if it is read-only, since
    #   we will remove it / modify it
    if Metis::File.exists?(new_file_name, new_bucket, new_folder)
      dest_file = get_file_obj_from_path(revision.dest)
      raise Etna::Forbidden, "Destination file #{dest_file.file_name} is read-only" if dest_file.read_only?
    end

    raise Etna::Forbidden, "Cannot copy over existing folder #{revision.dest}" if  Metis::Folder.exists?(new_file_name, new_bucket, new_folder)
  end
end
