class FilesController < Metis::Controller
  def list
    bucket = require_bucket

    folder = Metis::Folder.from_path(bucket, @params[:folder_path]).last

    if @params[:folder_path]
      raise Etna::Forbidden, 'No such folder!' unless folder
    end

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

    # any valid file_name is a valid folder name
    raise Etna::BadRequest, 'Invalid folder name' unless Metis::File.valid_file_path?(@params[:folder_path])

    # require the parent folder (if any) to exist
    parent_folder = nil
    parent_folder_path, folder_name = Metis::File.path_parts(@params[:folder_path])

    if parent_folder_path
      parent_folder = Metis::Folder.from_path(bucket, parent_folder_path).last

      raise Etna::BadRequest, 'Invalid parent folder' unless parent_folder
    end

    # is there a previous folder here?
    folder = Metis::Folder.where(
      project_name: @params[:project_name],
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
