require_relative '../../folder_rename_revision'

class FolderController < Metis::Controller
  def list
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    success_json(list_folder_contents(bucket, folder))
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
    
    folder.update(updated_at: Time.now)
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

    folders = mkdir_p(bucket, @params[:folder_path], @params[:project_name], Metis::File.author(@user))
    success_json(folders: [ folders.last.to_hash ])
  end

  def remove
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::BadRequest, 'Folder is read-only' if folder.read_only?

    raise Etna::Forbidden, 'Parent folder is read-only' if folder.folder&.read_only?

    raise Etna::BadRequest, 'Folder is not empty' unless folder.can_remove?

    response = { folders: [ folder.to_hash ] }

    # actually remove the folder
    folder.remove!

    return success_json(response)
  end

  def protect
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    folder.protect!

    success_json(folders: [ folder.to_hash ])
  end

  def unprotect
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::BadRequest, 'Folder is not protected' unless folder.read_only?

    folder.unprotect!

    success_json(folders: [ folder.to_hash ])
  end

  def rename
    require_param(:new_folder_path)

    source_bucket = require_bucket

    dest_bucket = source_bucket
    if @params[:new_bucket_name]
      dest_bucket = require_bucket(@params[:new_bucket_name])
    end

    revision = Metis::FolderRenameRevision.new({
      source: Metis::Path.path_from_parts(
        @params[:project_name],
        source_bucket.name,
        @params[:folder_path]
      ),
      dest: Metis::Path.path_from_parts(
        @params[:project_name],
        dest_bucket.name,
        @params[:new_folder_path]
      ),
      user: @user
    })

    revision.set_bucket(revision.source, [source_bucket, dest_bucket])
    revision.set_bucket(revision.dest, [source_bucket, dest_bucket])

    # we need the dest parent folder here, so we call File#path_parts
    source_folder_path, _ = Metis::File.path_parts(@params[:folder_path])
    dest_folder_path, _ = Metis::File.path_parts(@params[:new_folder_path])

    revision.set_folder(
      revision.source,
      [require_folder(source_bucket, source_folder_path)].compact)
    revision.set_folder(
      revision.dest,
      [require_folder(dest_bucket, dest_folder_path)].compact) if dest_folder_path

    revision.validate

    return failure(422, errors: revision.errors) unless revision.valid?

    return success_json(folders: [ revision.revise!.to_hash ])
  end

  protected

  def mkdir_p(bucket, folder_path, project_name, author)
    existing_folders = Metis::Folder.from_path(bucket, folder_path, allow_partial_match: true)
    folder_names = folder_path.split(%r!/!)

    folder_names.inject([]) do |parents, folder_name|
      existing = existing_folders.shift
      next (parents << existing) unless existing.nil?


      if Metis::File.exists?(folder_name, bucket, parents.last)
        raise Etna::BadRequest, "Cannot overwrite existing file"
      end

      begin
        parents << Metis::Folder.create(
          folder_name: folder_name,
          bucket_id: bucket&.id,
          folder_id: parents.last&.id,
          project_name: project_name,
          author: author,
      )
      rescue  Sequel::UniqueConstraintViolation => e
        ## Can occur if two simult requests get to the create line after both reading no existing_folders.
        ## Because of the uniq index constraint, this will occur for one of the requests, while the other should succeed.
        ## In that case, fall back to querying the other transaction's entry.
        ## Note: find_or_create does not fix this, it still does not handle the unique constraint and simply
        ## queries or creates, which is not good enough for READ COMMITTED isolation where a read might not see a yet
        ## committed create.
        Metis.instance.logger.log_error(e)
        parents << Metis::Folder.find(bucket_id: bucket&.id, folder_id: parents.last&.id, folder_name: folder_name)
      end
    end
  end

  def list_folder_contents(bucket, folder)
    files = Metis::File.where(
      project_name: @params[:project_name],
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |file|
      file.to_hash(request: @request)
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
