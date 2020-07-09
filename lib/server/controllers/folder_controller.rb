class FolderController < Metis::Controller
  def list
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    files = Metis::File.where(
      project_name: @params[:project_name],
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |file|
      file.to_hash(@request)
    end

    folders = Metis::Folder.where(
      bucket: bucket,
      folder_id: folder ? folder.id : nil
    ).all.map do |fold|
      fold.to_hash
    end

    success_json(files: files, folders: folders)
  end

  def create
    bucket = require_bucket
    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:folder_path])

    folders = mkdir_p(bucket, @params[:folder_path], @params[:project_name], Metis::File.author(@user))
    success_json(folders: [ folders.first.to_hash ])
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
    bucket = require_bucket
    folder = Metis::Folder.from_path(bucket, @params[:folder_path]).last

    raise Etna::Error.new('Folder not found', 404) unless folder

    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:new_folder_path])

    new_parent_folder_path, new_folder_name = Metis::File.path_parts(@params[:new_folder_path])

    new_parent_folder = require_folder(bucket, new_parent_folder_path)

    raise Etna::Forbidden, 'Folder is read-only' if new_parent_folder && new_parent_folder.read_only?

    raise Etna::BadRequest, 'Cannot overwrite existing folder' if Metis::Folder.exists?(new_folder_name, bucket, new_parent_folder)

    raise Etna::BadRequest, 'Cannot overwrite existing file' if Metis::File.exists?(new_folder_name, bucket, new_parent_folder)

    folder.rename!(new_parent_folder, new_folder_name)

    success_json(folders: [ folder.to_hash ])
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
        parents << Metis::Folder.find(bucket_id: bucket&.id, folder_id: parents.last&.id, folder_name: folder_name)
      end
    end
  end
end
