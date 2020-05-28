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

    revision = Metis::CopyRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        @params[:bucket_name],
        @params[:file_path]
      ),
      dest: Metis::Path.path_from_parts(
        @params[:project_name],
        @params[:new_bucket_name] ? @params[:new_bucket_name] : @params[:bucket_name],
        @params[:new_file_path]
      ),
      user: @user
    })

    revision.validate

    raise Etna::BadRequest, revision.errors unless revision.valid?

    success_json(files: [ revision.revise!.to_hash(@request) ])
  end

  def bulk_copy
    binding.pry
    # We want to support bulk copy operations
    #   from Magma.
    # First, check that the user has access to all of
    #   the revision buckets and files, and that
    #   all revisions are valid.
    # If all revisions are valid, execute them.
    require_param(:revisions)

    revisions = @params[:revisions].
      map {|rev| Metis::CopyRevision.new(rev.merge({
        user: @user
      })) }

    raise Etna::BadRequest, 'At least one revision required' unless revisions.length > 0

    # In bulk copy mode, we fetch the buckets and folders in bulk
    #   and set them on the revision objects, to try and minimize database hits
    hmac = @request.env['etna.hmac']

    user_authorized_buckets = Metis::Bucket.where(
      project_name: @params[:project_name],
      owner: ['metis', hmac.id.to_s],
      name: revisions.map(&:bucket_names).flatten.uniq
    ).all.select{|b| b.allowed?(@user, hmac)}

    revision_bucket_folders = {}

    revisions.each do |rev|
      rev.set_bucket(rev.source, user_authorized_buckets)
      rev.set_bucket(rev.dest, user_authorized_buckets)
    end

    user_authorized_buckets.each do |bucket|
      # Bulk-fetch all the unique folders that are in the
      #   revisions.
      bucket_folder_paths = revisions.
        map(&:paths).flatten.
        select{|p| p.bucket_name == bucket.name}.
        map(&:folder_path).flatten.unique

      bucket_folders = bucket_folder_paths.map {
        |folder_path| Metis::Folder.from_path(folder_path)
      }

      revisions.each do |rev|
        rev.set_folder(rev.source, bucket_folders)
        rev.set_folder(rev.dest, bucket_folders)
      end
    end

    revisions.map { |rev| rev.validate }

    errors = revisions.map(&:errors).flatten

    raise Etna::BadRequest, errors unless errors.length == 0

    # If we've gotten here, every revision looks good and we can execute them!
    success_json(
      files: revisions.map(&:revise!).
        map {|new_file| new_file.to_hash(@request)}
    )
  end
end
