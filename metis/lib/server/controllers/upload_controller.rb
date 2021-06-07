class UploadController < Metis::Controller
  def authorize
    require_params(:project_name, :bucket_name, :file_path)
    bucket = require_bucket

    raise Etna::BadRequest, 'Invalid path' unless Metis::File.valid_file_path?(@params[:file_path])

    folder_path, file_name = Metis::File.path_parts(@params[:file_path])

    folder = require_folder(bucket, folder_path)

    raise Etna::Forbidden, 'Folder is read-only' if folder&.read_only?

    raise Etna::Forbidden, 'Cannot overwrite existing folder' if Metis::Folder.exists?(file_name, bucket, folder)

    file = Metis::File.from_folder(bucket, folder, file_name)

    raise Etna::Forbidden, 'File cannot be overwritten' if file && file.read_only?

    # Make a MAC url
    url = Metis::File.upload_url(
      @request,
      @params[:project_name],
      @params[:bucket_name],
      @params[:file_path],
      @user
    )

    success_json(url: url)
  end

  UPLOAD_ACTIONS=[ :start, :blob, :cancel, :reset ]

  # this endpoint handles multiple possible actions, allowing us to authorize
  # one path /upload and support several upload operations
  def upload
    require_param(:action)

    action = @params[:action].to_sym

    raise Etna::BadRequest, 'Incorrect upload action' unless UPLOAD_ACTIONS.include?(action)

    send :"upload_#{action}"
  end

  private


  # create a metadata entry in the database and also a file on
  # the file system with 0 bytes.
  def upload_start
    require_params(:file_size, :next_blob_size, :next_blob_hash)
    bucket = require_bucket

    hmac = @request.env['etna.hmac']
    user = user_by_hmac(hmac)

    upload = Metis::Upload.fetch(
      project_name: @params[:project_name],
      file_name: @params[:file_path],
      bucket: bucket,
      metis_uid: metis_uid,
      user: user,
    )

    requires_reset = @params[:reset] || @params[:file_size] != upload.file_size

    # the upload has been started already, report the current
    # position
    if upload.current_byte_position > 0 && !requires_reset
      return success_json(upload)
    end

    upload_update = {
        file_size: @params[:file_size].to_i,
        next_blob_size: @params[:next_blob_size],
        next_blob_hash: @params[:next_blob_hash],
        current_byte_position: 0,
    }

    if requires_reset
      upload_update[:author] = Metis::File.author(user)
      upload.delete_partial!
    end

    upload.update(upload_update)

    # Send upload initiated
    success_json(upload)
  end

  # Upload a chunk of the file.
  def upload_blob
    require_params(:blob_data, :current_byte_position, :next_blob_size, :next_blob_hash)
    bucket = require_bucket

    upload = Metis::Upload.where(
      project_name: @params[:project_name],
      file_name: @params[:file_path],
      bucket: bucket,
      metis_uid: metis_uid
    ).first

    raise Etna::BadRequest, 'Upload has not been started' unless upload

    raise Etna::BadRequest, 'Wrong byte position' unless @params[:current_byte_position].to_i == upload.current_byte_position

    blob = Metis::Blob.new(@params[:blob_data])

    raise Etna::BadRequest, 'Blob integrity failed' unless blob.continues?(upload)

    upload.append_blob(
      blob,
      @params[:next_blob_size],
      @params[:next_blob_hash]
    )

    return complete_upload(upload) if upload.complete?

    return success_json(upload)
  end

  private

  def user_by_hmac(hmac)
    return Etna::User.new(
        email: (hmac.headers[:email] || hmac.id).to_s,
        name: (hmac.headers[:name] || hmac.id).to_s) if hmac and hmac.valid?
    @user
  end

  def complete_upload(upload)
    folder_path, file_name = Metis::File.path_parts(upload.file_name)
    folder = require_folder(upload.bucket, folder_path)

    if folder && folder.read_only?
      raise Etna::Forbidden, 'Folder is read-only'
    end

    file = Metis::File.from_folder(upload.bucket, folder, file_name)

    if file && file.read_only?
      raise Etna::Forbidden, 'Cannot overwrite existing file'
    end

    if Metis::Folder.exists?(file_name, upload.bucket, folder)
      raise Etna::Forbidden, 'Cannot overwrite existing folder'
    end

    new_file = upload.finish!

    # we will embed the new file hash inside the
    # upload hash
    upload_hash = upload.to_hash
    upload_hash[:file] = new_file.to_hash(request: @request)

    return success_json(upload_hash)
  ensure
    upload.delete_with_partial!
  end

  public

  def upload_cancel
    bucket = require_bucket

    upload = Metis::Upload.where(
      project_name: @params[:project_name],
      file_name: @params[:file_path],
      bucket: bucket,
      metis_uid: metis_uid
    ).first

    raise Etna::BadRequest, 'Upload has not been started' unless upload

    response = success_json(upload: upload.to_hash)
    # axe the upload data and record
    upload.delete_with_partial!

    return response
  end

  private

  def metis_uid
    id = @request.cookies[Metis.instance.config(:metis_uid_name)] || @params[:metis_uid]

    if id.nil?
      raise Etna::BadRequest("metis_uid not set, did you forget to configure metis_uid to #{Metis.instance.config(:metis_uid_name)}?")
    end

    id
  end
end
