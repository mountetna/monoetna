require_relative '../../folder_rename_revision'
require_relative '../../folder_copy_revision'

class FolderController < Metis::Controller
  def list
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    success_json(list_folder_contents(bucket, folder))
  end

  def size
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    success_json(
      size: folder ? folder.size : bucket.size
    )
  end

  def list_by_id
    bucket = require_bucket
    folder = require_folder_by_id(bucket, @params[:folder_id])

    success_json(list_folder_contents(bucket, folder))
  end

  def touch
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::BadRequest, "Invalid folder: #{@params[:folder_path]}" unless folder

    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    folder.update(updated_at: Time.now, author: Metis::File.author(@user))
    folder.refresh

    success_json(folders: [ folder.to_hash ])
  end

  def list_all_folders
    bucket = require_bucket

    limit = @params.has_key?(:limit) ? @params[:limit].to_i : nil
    offset = @params.has_key?(:offset) ? @params[:offset].to_i : nil

    raise Etna::BadRequest, "Invalid offset" if offset&.negative?
    raise Etna::BadRequest, "Invalid limit" if limit&.negative?

    folders = Metis::Folder.where(
      bucket: bucket
    ).all

    limit = limit ? limit : folders.length
    offset = offset ? offset : 0

    folder_hashes = folder_hashes_with_calculated_paths(
      all_folders: folders,
      offset: offset,
      limit: limit,
      bucket: bucket)

    success_json(
      folders: folder_hashes)
  end

  def create
    bucket = require_bucket
    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:folder_path])

    folders = Metis::Folder.mkdir_p(bucket, @params[:folder_path], @params[:project_name], Metis::File.author(@user))

    event_log(
      event: 'create_folder',
      message: "created folders in bucket #{@params[:bucket_name]}",
      payload: {
        folders: [ @params[:folder_path] ]
      },
      consolidate: true
    )

    success_json(folders: [ folders.last.to_hash ])
  end

  def remove
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::BadRequest, 'Folder is read-only' if folder.read_only?

    raise Etna::Forbidden, 'Parent folder is read-only' if folder.folder&.read_only?

    raise Etna::BadRequest, 'Folder is not empty' unless @params[:recursive] || folder.can_remove?

    response = { folders: [ folder.to_hash ] }

    # remove contents if necessary
    if @params[:recursive]
      folder.remove_contents!
    end

    # actually remove the folder
    folder.remove!

    event_log(
      event: 'remove_folder',
      message: "removed folders in bucket #{@params[:bucket_name]}",
      payload: {
        folders: [ @params[:folder_path] ]
      },
      consolidate: true
    )

    return success_json(response)
  end

  def protect
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    folder.protect!

    event_log(
      event: 'protect_folder',
      message: "protected folders in bucket #{@params[:bucket_name]}",
      payload: {
        folders: [ @params[:folder_path] ]
      },
      consolidate: true
    )

    success_json(folders: [ folder.to_hash ])
  end

  def unprotect
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::BadRequest, 'Folder is not protected' unless folder.read_only?

    folder.unprotect!

    event_log(
      event: 'unprotect_folder',
      message: "unprotected folders in bucket #{@params[:bucket_name]}",
      payload: {
        folders: [ @params[:folder_path] ]
      },
      consolidate: true
    )

    success_json(folders: [ folder.to_hash ])
  end

  def rename
    require_param(:new_folder_path)

    source_bucket = require_bucket

    dest_bucket = source_bucket
    if @params[:new_bucket_name]
      dest_bucket = require_bucket(@params[:new_bucket_name])
    end

    new_folder_path = Metis::Path.remove_double_slashes(@params[:new_folder_path])

    revision = Metis::FolderRenameRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        source_bucket.name,
        @params[:folder_path]
      ),
      dest: Metis::Path.path_from_parts(
        @params[:project_name],
        dest_bucket.name,
        new_folder_path
      ),
      user: @user
    })

    revision.set_bucket(revision.source, [source_bucket, dest_bucket])
    revision.set_bucket(revision.dest, [source_bucket, dest_bucket])

    # we need the dest parent folder here, so we call File#path_parts
    source_folder_path, _ = Metis::File.path_parts(@params[:folder_path])
    dest_folder_path, _ = Metis::File.path_parts(new_folder_path)

    revision.set_folder(
      revision.source,
      [require_folder(source_bucket, source_folder_path)].compact)
    revision.set_folder(
      revision.dest,
      [require_folder(dest_bucket, dest_folder_path)].compact) if dest_folder_path

    revision.validate

    return failure(422, errors: revision.errors) unless revision.valid?

    event_log(
      event: 'rename_folder',
      message: "renamed folders in bucket #{@params[:bucket_name]}",
      payload: {
        folders: [ @params[:folder_path] ]
      },
      consolidate: true
    )

    return success_json(folders: [ revision.revise!.to_hash ])
  end

  def copy
    require_param(:new_parent_folder)

    source_bucket = require_bucket

    dest_bucket = source_bucket
    if @params[:new_bucket_name]
      dest_bucket = require_bucket(@params[:new_bucket_name])
    end

    new_parent_folder_path = Metis::Path.remove_double_slashes(@params[:new_parent_folder])

    revision = Metis::FolderCopyRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        source_bucket.name,
        @params[:folder_path]
      ),
      dest: Metis::Path.path_from_parts(
        @params[:project_name],
        dest_bucket.name,
        new_parent_folder_path
      ),
      user: @user
    })

    revision.set_bucket(revision.source, [source_bucket, dest_bucket])
    revision.set_bucket(revision.dest, [source_bucket, dest_bucket])

    # Don't' use set_folder because that will look for the parent folders as
    #   if the path passed in is for a file ...
    #   but in this case, the user is passing in folder paths already.
    revision.source.folder = require_folder(source_bucket, @params[:folder_path])
    revision.dest.folder = require_folder(dest_bucket, new_parent_folder_path) unless new_parent_folder_path.empty?

    revision.validate

    return failure(422, errors: revision.errors) unless revision.valid?

    event_log(
      event: 'copy_folder',
      message: "copied folders in bucket #{@params[:bucket_name]}",
      payload: {
        folders: [ @params[:folder_path] ]
      },
      consolidate: true
    )

    return success_json(folders: [ revision.revise!.to_hash ])
  end

  protected

  def list_folder_contents(bucket, folder)
    files = Metis::File.where(
      project_name: @params[:project_name],
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |file|
      file.to_hash(user: @user)
    end

    folders = Metis::Folder.where(
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |fold|
      fold.to_hash
    end

    {
      files: files,
      folders: folders
    }
  end
end
