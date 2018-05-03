class UploadController < Metis::Controller
  def authorize
    require_params(:project_name, :file_name)

    raise Etna::BadRequest, 'The file name is illegal.' unless Metis::File.valid_file_name?(@params[:file_name])

    bucket_name = @params[:bucket_name] || 'files'

    bucket = Metis::Bucket.find(name: bucket_name, project_name: @params[:project_name])

    raise Etna::BadRequest, 'No such bucket!' unless bucket

    raise Etna::Forbidden, 'Inaccessible bucket.' unless bucket.allowed?(@user)

    file = Metis::File.find(
      project_name: @params[:project_name],
      file_name: @params[:file_name],
      bucket: bucket
    )

    raise Etna::Forbidden, 'File cannot be overwritten.' if file && file.read_only?

    # Create the upload
    upload = Metis::Upload.create(
      file_name: @params[:file_name],
      project_name: @params[:project_name],
      author: Metis::File.author(@user),
      bucket: bucket,
      metis_uid: @request.cookies[Metis.instance.config(:metis_uid_name)],
      file_size: 0,
      current_byte_position: 0,
      next_blob_size: -1,
      next_blob_hash: ''
    )

    # Make a MAC url
    url = Metis::File.upload_url(
      @request,
      @params[:project_name],
      bucket_name,
      @params[:file_name]
    )

    success(url)
  end

  UPLOAD_ACTIONS=[ :start, :blob, :cancel, :reset ]

  # this endpoint handles multiple possible actions, allowing us to authorize
  # one path /upload and support several upload operations
  def upload
    require_params(:project_name, :file_name, :action)

    action = @params[:action].to_sym

    raise Etna::BadRequest, 'Incorrect upload action' unless UPLOAD_ACTIONS.include?(action)

    send :"upload_#{action}"
  end

  private

  # create a metadata entry in the database and also a file on
  # the file system with 0 bytes.
  def upload_start
    require_params(:file_size, :next_blob_size, :next_blob_hash)

    # get the current bucket
    bucket = Metis::Bucket.find(name: @params[:bucket_name])

    upload = Metis::Upload.where(
      project_name: @params[:project_name],
      file_name: @params[:file_name],
      bucket: bucket,
      metis_uid: @request.cookies[Metis.instance.config(:metis_uid_name)]
    ).first

    raise Etna::BadRequest, 'No matching upload!' unless upload

    # the upload has been started already, report the current
    # position
    if upload.current_byte_position > 0
      return success_json(upload)
    end

    upload.update(
      file_size: @params[:file_size].to_i,
      next_blob_size: @params[:next_blob_size],
      next_blob_hash: @params[:next_blob_hash]
    )

    # Send upload initiated
    success_json(upload)
  end

  # Upload a chunk of the file.
  def upload_blob
    require_params(:blob_data, :next_blob_size, :next_blob_hash)

    bucket = Metis::Bucket.find(name: @params[:bucket_name])

    upload = Metis::Upload.where(
      project_name: @params[:project_name],
      file_name: @params[:file_name],
      bucket: bucket,
      metis_uid: @request.cookies[Metis.instance.config(:metis_uid_name)],
    ).first

    raise Etna::BadRequest, 'Upload has not been started!' unless upload

    blob_path = @params[:blob_data][:tempfile].path

    raise Etna::BadRequest, 'Blob integrity failed' unless upload.blob_valid?(blob_path)

    upload.append_blob(blob_path)

    upload.update(
      next_blob_size: @params[:next_blob_size],
      next_blob_hash: @params[:next_blob_hash]
    )

    if upload.complete?
      if upload.can_place?
        upload.finish!

        response = success_json(upload)

        upload.delete

        return response
      else
        upload.delete_partial!
        upload.delete

        raise Etna::Forbidden, 'Cannot overwrite existing file!'
      end
    end

    return success_json(upload)
  end

  def upload_cancel
    bucket = Metis::Bucket.find(name: @params[:bucket_name])

    upload = Metis::Upload.where(
      project_name: @params[:project_name],
      file_name: @params[:file_name],
      bucket: bucket,
      metis_uid: @request.cookies[Metis.instance.config(:metis_uid_name)],
    ).first

    raise Etna::BadRequest, 'Upload has not been started!' unless upload

    # axe the upload data and record
    upload.delete_partial!
    upload.delete

    return success('deleted')
  end
end
