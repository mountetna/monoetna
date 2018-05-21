class FilesController < Metis::Controller
  def list
    bucket = Metis::Bucket.where(project_name: @params[:project_name], name: @params[:bucket_name]).first

    raise Etna::BadRequest, 'Invalid bucket!' unless bucket && bucket.allowed?(@user)

    folder = nil

    if @params[:folder_name]
      folder = Metis::File.where(
        project_name: @params[:project_name],
        file_name: @params[:folder_name],
        is_folder: true
      ).first
      raise Etna::Forbidden, 'No such folder!' unless folder
    end

    files = Metis::File.where(
      project_name: @params[:project_name],
      folder_id: folder ? folder.id : nil
    ).all.map do |file|
      file.to_hash(@request)
    end

    success({ files: files }.to_json, 'application/json')
  end

  def create_folder
    bucket = Metis::Bucket.where(project_name: @params[:project_name], name: @params[:bucket_name]).first

    raise Etna::BadRequest, 'Invalid bucket!' unless bucket && bucket.allowed?(@user)

    # any valid file_name is a valid folder name
    raise Etna::BadRequest, 'Invalid folder name' unless Metis::File.valid_file_name?(@params[:folder_name])

    # require the parent folder (if any) to exist
    parent_folder = nil
    parent_folder_name = Metis::File.foldername(@params[:folder_name])

    if parent_folder_name
      parent_folder = Metis::File.find(
        project_name: @params[:project_name],
        bucket: bucket,
        file_name: parent_folder_name
      )

      raise Etna::BadRequest, 'Invalid parent folder' unless parent_folder && parent_folder.folder?
    end

    # is there a previous file here?
    file = Metis::File.where(project_name: @params[:project_name], bucket: bucket, file_name: @params[:folder_name]).first
    raise Etna::BadRequest, "#{file.folder? ? 'Folder' : 'File'} exists" if file

    # create the folder
    file = Metis::File.create(
      project_name: @params[:project_name],
      file_name: @params[:folder_name],
      author: Metis::File.author(@user),
      bucket: bucket,
      folder: parent_folder,
      is_folder: true
    )

    success({ folders: [ file.to_hash(@request) ] })
  end
end
