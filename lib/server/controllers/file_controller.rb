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
    bucket = require_bucket
    file = Metis::File.from_path(bucket, @params[:file_path])

    raise Etna::Error.new('File not found', 404) unless file&.has_data?

    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:new_file_path])

    raise Etna::Forbidden, "Cannot edit project #{@params[:project_name]}" unless @user.can_edit?(@params[:project_name])

    new_folder_path, new_file_name = Metis::File.path_parts(@params[:new_file_path])

    new_folder = require_folder(bucket, new_folder_path)

    raise Etna::Forbidden, 'Folder is read-only' if new_folder&.read_only?

    raise Etna::Forbidden, 'Cannot copy over existing file' if Metis::File.exists?(new_file_name, bucket, new_folder)

    raise Etna::Forbidden, 'Cannot copy over existing folder' if  Metis::Folder.exists?(new_file_name, bucket, new_folder)

    new_file = Metis::File.create(
      project_name: @params[:project_name],
      file_name: new_file_name,
      folder_id: new_folder&.id,
      bucket: bucket,
      author: Metis::File.author(@user),
      data_block: file.data_block
    )
    success_json(files: [ new_file.to_hash(@request) ])
  end
end
