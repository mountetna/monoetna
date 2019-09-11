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

    parent_folder_path, folder_name = Metis::File.path_parts(@params[:folder_path])
    parent_folder = require_folder(bucket, parent_folder_path)

    # is there a previous folder here?
    folder = Metis::Folder.where(
      bucket: bucket,
      folder_id: parent_folder ? parent_folder.id : nil,
      folder_name: folder_name
    ).first
    raise Etna::BadRequest, "Folder exists" if folder

    # create the folder
    folder = Metis::Folder.create(
      project_name: @params[:project_name],
      folder_name: folder_name,
      author: Metis::File.author(@user),
      bucket: bucket,
      folder: parent_folder
    )
    folder.create_actual_folder!

    success_json(folders: [ folder.to_hash ])
  end

  def remove
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    raise Etna::BadRequest, 'Folder is read-only' if folder.read_only?

    raise Etna::BadRequest, 'Folder is not empty' unless folder.can_remove?

    response = success_json(folders: [ folder.to_hash ])

    # actually remove the folder
    folder.remove!

    return response
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

    raise Etna::Error.new('Folder not found', 404) unless folder && folder.has_directory?

    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:new_folder_path])

    new_parent_folder_path, new_folder_name = Metis::File.path_parts(@params[:new_folder_path])

    new_parent_folder = require_folder(bucket, new_parent_folder_path)

    raise Etna::Forbidden, 'Folder is read-only' if new_parent_folder && new_parent_folder.read_only?

    existing_new_folder = Metis::Folder.where(
      bucket: bucket,
      folder_id: new_parent_folder ? new_parent_folder.id : nil,
      folder_name: new_folder_name).first

    raise Etna::Forbidden, 'Cannot overwrite existing folder' if existing_new_folder

    folder.rename!(new_parent_folder, new_folder_name)

    success_json(folders: [ folder.to_hash ])
  end
end
