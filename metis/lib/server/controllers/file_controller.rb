require_relative '../../copy_revision'
require_relative '../../file_rename_revision'

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

    success_json(files: [ file.to_hash(request: @request) ])
  end

  def unprotect
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'File is not protected' unless file.read_only?

    file.unprotect!

    success_json(files: [ file.to_hash(request: @request) ])
  end

  def touch
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::Forbidden, 'File is read only' if file.read_only?
    
    file.update(updated_at: Time.now)
    file.refresh

    success_json(files: [ file.to_hash ])
  end

  def rename
    require_param(:new_file_path)

    source_bucket = require_bucket

    dest_bucket = source_bucket
    if @params[:new_bucket_name]
      dest_bucket = require_bucket(@params[:new_bucket_name])
    end

    revision = Metis::FileRenameRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        source_bucket.name,
        @params[:file_path]
      ),
      dest: Metis::Path.path_from_parts(
        @params[:project_name],
        dest_bucket.name,
        @params[:new_file_path]
      ),
      user: @user
    })

    revision.set_bucket(revision.source, [source_bucket, dest_bucket])
    revision.set_bucket(revision.dest, [source_bucket, dest_bucket])

    # we need the dest parent folder here, so we call File#path_parts
    source_folder_path, _ = Metis::File.path_parts(@params[:file_path])
    dest_folder_path, _ = Metis::File.path_parts(@params[:new_file_path])

    revision.set_folder(
      revision.source,
      [require_folder(source_bucket, source_folder_path)].compact)
    revision.set_folder(
      revision.dest,
      [require_folder(dest_bucket, dest_folder_path)].compact) if dest_folder_path

    revision.validate

    return failure(422, errors: revision.errors) unless revision.valid?

    return success_json(files: [ revision.revise!.to_hash ])
  end

  def copy
    require_param(:new_file_path)

    source_bucket = require_bucket

    dest, dest_bucket = construct_dest_info

    revision = Metis::CopyRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        source_bucket.name,
        @params[:file_path]
      ),
      dest: dest.path,
      user: @user
    })

    revision.set_bucket(revision.source, [source_bucket, dest_bucket])
    revision.set_bucket(revision.dest, [source_bucket, dest_bucket])

    source_folder_path, _ = Metis::File.path_parts(@params[:file_path])

    revision.set_folder(
      revision.source,
      [require_folder(source_bucket, source_folder_path)]) if source_folder_path
    revision.set_folder(
      revision.dest,
      [require_folder(dest_bucket, dest.folder_path)]) if dest.folder_path

    revision.validate

    return failure(422, errors: revision.errors) unless revision.valid?

    return success_json(files: [ revision.revise!.to_hash(request: @request) ])
  end

  def bulk_copy
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
    user_authorized_buckets = Metis::Bucket.where(
      owner: ['metis', @hmac&.id.to_s].reject {|o| o.empty?},
      name: revisions.map(&:bucket_names).flatten.uniq
    ).all.select{|b| b.allowed?(@user, @hmac)}

    revision_bucket_folders = {}

    revisions.each do |rev|
      rev.set_bucket(rev.source, user_authorized_buckets)
      rev.set_bucket(rev.dest, user_authorized_buckets)
    end

    user_authorized_buckets.each do |bucket|
      # Bulk-fetch all the unique folders that are in the
      #   revisions.
      bucket_folder_paths = revisions.
        map(&:mpaths).flatten.
        select{|p| p.bucket_matches?(bucket)}.
        map(&:folder_path).flatten.compact.uniq

      bucket_folders = bucket_folder_paths.map {
        |folder_path| Metis::Folder.from_path(bucket, folder_path).last
      }.flatten.compact

      revisions.each do |rev|
        if rev.source.mpath.bucket_matches?(bucket)
          rev.set_folder(rev.source, bucket_folders)
        end
        if rev.dest.mpath.bucket_matches?(bucket)
          rev.set_folder(rev.dest, bucket_folders)
        end
      end
    end

    revisions.map { |rev| rev.validate }

    errors = revisions.select{|r| !r.valid?}.map(&:to_hash)

    return failure(422, errors: errors) unless errors.length == 0

    # If we've gotten here, every revision looks good and we can execute them!
    return success_json(
      files: revisions.map(&:revise!).
        map {|new_file| new_file.to_hash(request: @request)}
    )
  end

  private

  def construct_dest_info
    if Metis::Path.filepath_match.match(@params[:new_file_path])
      dest = Metis::Path.new(@params[:new_file_path])
      dest_bucket = require_bucket(dest.bucket_name, dest.project_name)
      return dest, dest_bucket
    end

    dest_bucket = require_bucket
    if @params[:new_bucket_name]
      dest_bucket = require_bucket(@params[:new_bucket_name])
    end

    dest_path = Metis::Path.path_from_parts(
      @params[:project_name],
      dest_bucket.name,
      @params[:new_file_path]
    )

    return Metis::Path.new(dest_path), dest_bucket
  end
end
