require_relative '../../copy_revision'
require_relative '../../file_rename_revision'

class FileController < Metis::Controller
  def remove
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::Forbidden, 'File is restricted' if file.restrict_user?(@user)

    raise Etna::Forbidden, 'File is read-only' if file.read_only?

    raise Etna::Forbidden, 'Folder is read-only' if file.folder&.read_only?

    response = { files: [ file.to_hash ] }

    # Capture file path and datablock before deletion
    file_path = file.file_path
    datablock = file.data_block

    file.remove!

    # Log unlink event AFTER successful deletion
    Metis::DataBlockLedger.log_unlink(file, datablock, @user, file_path: file_path)

    event_log(
      event: 'remove_file',
      message: "removed files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ @params[:file_path] ]
      },
      consolidate: true
    )
    return success_json(response)
  end

  def restrict
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'File is already restricted' if file.restricted?

    file.restrict!

    event_log(
      event: 'restrict_file',
      message: "restricted files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ @params[:file_path] ]
      },
      consolidate: true
    )

    success_json(files: [ file.to_hash(user: @user) ])
  end

  def unrestrict
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'File is not restricted' unless file.restricted?

    file.unrestrict!

    event_log(
      event: 'unrestrict_file',
      message: "unrestricted files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ @params[:file_path] ]
      },
      consolidate: true
    )

    success_json(files: [ file.to_hash(user: @user) ])
  end

  def protect
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::Forbidden, 'File is restricted' if file.restrict_user?(@user)

    raise Etna::BadRequest, 'File is already read-only' if file.read_only?

    file.protect!

    event_log(
      event: 'protect_file',
      message: "protected files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ @params[:file_path] ]
      },
      consolidate: true
    )

    success_json(files: [ file.to_hash(user: @user) ])
  end

  def unprotect
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::Forbidden, 'File is restricted' if file.restrict_user?(@user)

    raise Etna::BadRequest, 'File is not protected' unless file.read_only?

    file.unprotect!

    event_log(
      event: 'unprotect_file',
      message: "unprotected files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ @params[:file_path] ]
      },
      consolidate: true
    )

    success_json(files: [ file.to_hash(user: @user) ])
  end

  def touch_files
    require_param(:file_paths)
    bucket = require_bucket

    files = @params[:file_paths].uniq.map do |file_path|
      [ file_path, Metis::File.from_path(bucket, file_path) ]
    end

    files, missing_files = files.partition do |file_path, file|
      file&.has_data?
    end.map(&:to_h)

    return failure(404, errors: missing_files.map do |file_path, file|
      "File not found #{file_path}"
    end) unless missing_files.empty?

    files, restricted_files = files.partition do |file_path, file|
      !file.restrict_user?(@user)
    end.map(&:to_h)

    return failure(403, errors: restricted_files.map do |file_path, file|
      "File is restricted #{file_path}"
    end) unless restricted_files.empty?

    files, read_only_files = files.partition do |file_path, file|
      !file.read_only?
    end.map(&:to_h)

    return failure(403, errors: read_only_files.map do |file_path, file|
      "File is read-only #{file_path}"
    end) unless read_only_files.empty?

    Metis::File.where(
      id: files.values.map(&:id)
    ).update(
      updated_at: Time.now,
      author: Metis::File.author(@user)
    )

    files = Metis::File.where(
      id: files.values.map(&:id)
    ).all

    success_json(files: files.map(&:to_hash))
  end

  def touch
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::NotFound, 'File not found' unless file&.has_data?

    raise Etna::Forbidden, 'File is restricted' if file.restrict_user?(@user)

    raise Etna::Forbidden, 'File is read only' if file.read_only?

    file.update(updated_at: Time.now, author: Metis::File.author(@user))
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

    event_log(
      event: 'rename_file',
      message: "renamed files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ { old_file_path: @params[:file_path], new_file_path: @params[:new_file_path] } ]
      },
      consolidate: true
    )

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

    event_log(
      event: 'copy_file',
      message: "copied files in bucket #{@params[:bucket_name]}",
      payload: {
        files: [ { old_file_path: @params[:file_path], new_file_path: @params[:new_file_path] } ]
      },
      consolidate: true
    )

    return success_json(files: [ revision.revise!.to_hash(user: @user) ])
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

    read_only_buckets = reserved_buckets

    revisions.each do |rev|
      rev.set_bucket(rev.source, user_authorized_buckets + read_only_buckets)
      rev.set_bucket(rev.dest, user_authorized_buckets)
    end

    set_revision_folder(revisions, user_authorized_buckets, set_dest: true)
    set_revision_folder(revisions, read_only_buckets, set_dest: false)

    revisions.map { |rev| rev.validate }

    errors = revisions.select{|r| !r.valid?}.map(&:to_hash)

    return failure(422, errors: errors) unless errors.length == 0

    event_log(
      event: 'copy_file',
      message: "bulk-copied files in buckets #{@params[:bucket_name]}",
      payload: {
        files: revisions.map do |revision|
          { source: revision.source.mpath.path, dest: revision.dest.mpath.path }
        end
      },
      consolidate: true
    )
    # If we've gotten here, every revision looks good and we can execute them!
    return success_json(
      files: revisions.map(&:revise!).
        map {|new_file| new_file.to_hash(user: @user)}
    )
  end

  private

  def reserved_buckets
    Metis::Bucket.reserved_buckets_for_project(@params[:project_name])
  end

  def set_revision_folder(revisions, buckets, set_dest: false)
    buckets.each do |bucket|
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
        if set_dest && rev.dest.mpath.bucket_matches?(bucket)
          rev.set_folder(rev.dest, bucket_folders)
        end
      end
    end
  end

  def construct_dest_info
    new_file_path = Metis::Path.remove_double_slashes(@params[:new_file_path])

    if Metis::Path.filepath_match.match(new_file_path)
      dest = Metis::Path.new(new_file_path)
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
      new_file_path
    )

    return Metis::Path.new(dest_path), dest_bucket
  end
end
