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

    success_json(folders: [ folder.to_hash ])
  end
end
