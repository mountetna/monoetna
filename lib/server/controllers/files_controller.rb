class FilesController < Metis::Controller
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

  def create_folder
    bucket = require_bucket
    parent_folder_path, folder_name = parse_path(@params[:folder_path])
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

  def remove_folder
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::BadRequest, 'Cannot remove folder' unless folder.can_remove?

    response = success_json(folders: [ folder.to_hash ])

    folder.remove!

    return response
  end

  def protect_folder
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::Forbidden, 'Folder is read-only' if folder.read_only?

    folder.protect!

    success_json(folders: [ folder.to_hash ])
  end

  def unprotect_folder
    bucket = require_bucket
    folder = require_folder(bucket, @params[:folder_path])

    # remove the folder
    raise Etna::BadRequest, 'Folder is not protected' unless folder.read_only?

    folder.unprotect!

    success_json(folders: [ folder.to_hash ])
  end

  def remove_file
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    # remove the file
    raise Etna::BadRequest, 'Cannot remove file' unless file.can_remove?

    response = success_json(files: [ file.to_hash ])

    file.remove!

    return response
  end

  def protect_file
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    # remove the file
    raise Etna::Forbidden, 'File is read-only' if file.read_only?

    file.protect!

    success_json(files: [ file.to_hash ])
  end

  def unprotect_file
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file && file.has_data?

    # remove the file
    raise Etna::BadRequest, 'File is not protected' unless file.read_only?

    file.unprotect!

    success_json(files: [ file.to_hash ])
  end
end
