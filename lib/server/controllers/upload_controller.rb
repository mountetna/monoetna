class UploadController < Metis::Controller
  def authorize
    require_params(:project_name, :file_name)

    hmac = Etna::Hmac.new(
      Metis.instance,
      method: 'POST',
      host: @request.host,
      path: '/upload',
      expiration: (Time.now + Metis.instance.config(:upload_expiration)).iso8601,
      nonce: Metis.instance.sign.uid,
      id: :metis,
      headers: { project_name: @params[:project_name], file_name: @params[:file_name] }
    )

    url = URI::HTTPS.build(hmac.url_params)

    success(url)
  end

  UPLOAD_ACTIONS=[ :start, :blob ]

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

    # make an entry for the file if it does not exist
    file = Metis::File.find_or_create(
      project_name: @params[:project_name],
      file_name: @params[:file_name]
    ) do |f|
      f.original_name = @params[:file_name]
      f.uploader = ''
      f.size = 0
    end

    upload = Metis::Upload.where(
      file: file,
      metis_uid: @request.cookies[Metis.instance.config(:metis_uid_name)]
    ).first

    if upload
      return success(upload.to_json, 'application/json')
    end

    raise Etna::BadRequest, 'Upload in progress' if !file.uploads.empty?

    upload = Metis::Upload.create(
      file: file,
      status: 'initialized',
      metis_uid: @request.cookies[Metis.instance.config(:metis_uid_name)],
      file_size: @params[:file_size].to_i,
      current_byte_position: 0,
      current_blob_size: 0,
      next_blob_size: @params[:next_blob_size],
      next_blob_hash: @params[:next_blob_hash]
    )

    # Send upload initiated
    success(upload.to_json,'application/json')
  end

  # Upload a chunk of the file.
  def upload_blob
    require_params(:blob_data, :next_blob_size, :next_blob_hash)

    file = Metis::File.find(
      project_name: @params[:project_name],
      file_name: @params[:file_name]
    )

    raise Etna::BadRequest, 'Could not find file!' unless file

    upload = Metis::Upload.where(
      file: file,
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
      upload.finish!

      upload_json = upload.to_json

      upload.delete

      return success(upload_json, 'application/json')
    end

    return success(upload.to_json, 'application/json')
  end

  def pause_upload
    common_error_check
    @upload.update(status: 'paused')
    return send_upload_status
  end

  def cancel_upload
    common_error_check
    @upload.delete if @upload.respond_to?(:delete)
    @file.delete if @file.respond_to?(:delete) 
    remove_temp_file
    @params.delete('blob')
    @params[:status] = 'cancelled'
    return { success: true, request: @params }
  end

  def remove_file
    remove_file_error_check
    @upload.delete if @upload.respond_to?(:delete) 
    @file.delete if @file.respond_to?(:delete)
    remove_file_on_disk
    remove_temp_file
    @params['status'] = 'removed'
    return { :success=> true, :request=> @params }
  end

  def remove_failed
    remove_failed_error_check
    @upload.delete if @upload.respond_to?(:delete)
    @file.delete if @file.respond_to?(:delete)
    remove_temp_file
    @params['status'] = 'removed'
    return { :success=> true, :request=> @params }
  end

  def recover_upload
    recover_upload_error_check

    raise Etna::BadRequest, 'Missing param last_kib_hash!' unless @params[:last_kib_hash]

    @upload.update(
      status: 'paused',
      next_blob_size: 1024,
      next_blob_hash: @params[:last_kib_hash]
    )
    @params[:hashing_algorithm] = 'MD5'
    @params[:start_timestamp] = @file.start_upload
    @params[:token] = @user.token
    @params[:user_email] = @user.email
    @params[:hmac_signature] = generate_hmac

    @params.delete('first_kib_hash')
    @params.delete('last_kib_hash')
    @params['status'] = 'paused'
    return { :success=> true, :request=> @params }
  end

  private


  # Some common things we should check before we make an action.
  def common_error_check

    # Check the HMAC
    raise Etna::BadRequest, 'Invalid HMAC' unless hmac_valid?

    # Check that the file metadata exists and set global vars if they exist.
    file_metadata_check_and_set
  end

  # All the things we should check before we remove a file.
  def remove_file_error_check
    # Check that the user has edit permission on the project requested.
    user_edit_check

    get_file_metadata
    get_partial_metadata
  end

  def remove_failed_error_check

    # Check that the file metadata exists and set global vars if they exist.
    file_metadata_check_and_set

    # The file should NOT exist but the partial should.
    raise Etna::BadRequest, 'File exists' if file_exists?
    raise Etna::BadRequest, 'Partial does not exist' unless partial_exists?

    # Check that the user requesting the 'remove' is the same as the uploader.
    raise Etna::BadRequest, 'Deletion must be done by the uploader' unless @user.email == @file.upload_by
  end

  # All the items to check before we 'recover' an upload sequence.
  def recover_upload_error_check

    # Check that the correct parameters are present.
    raise Etna::BadRequest, 'Not authorized' unless has_auth_params?(@params)

    # Check that the user has edit permission on the project requested.
    user_edit_check

    # Check that this file system is in sync with the auth server. If there is a
    # group/project set in Janus there should be a corresponding directory
    # in Metis.
    raise Etna::BadRequest, 'Project dir missing' unless directory_exists?

    # Check that a full file does not exist. A partial should exist.
    raise Etna::BadRequest, "File exists" if file_exists?
    raise Etna::BadRequest, "Partial does not exist" unless partial_exists?

    # Check that the file metadata exists and set global vars if they exist.
    file_metadata_check_and_set

    # Check the first kilobyte hash.
    first_kilobyte_check
  end

  # Check that the file metadata exists and set global vars if they exist.
  # This method also sets the variables @file and @update. We need to find a
  # cleaner way to initialize these variables.
  def file_metadata_check_and_set

    get_file_metadata
    get_partial_metadata

    # The 'main' portion of the metadata should exist.
    raise_err(:BAD_REQ, 7, __method__) if !@file

    # The 'upload' portion of the metadata should exist.
    raise_err(:BAD_REQ, 7, __method__) if !@upload
  end

  # Check that the user has edit permission on the project requested.
  def user_edit_check

    if !@user.project_editor?(@params['project_name']) &&
       !@user.project_admin?(@params['project_name'])

      raise_err(:BAD_REQ, 8, __method__)
    end
  end

  def get_file_metadata
    @file = Metis::File[
      group_name: @params['group_name'],
      project_name: @params['project_name'],
      file_name: @params['file_name']
    ]
  end

  def get_partial_metadata

    if @file.respond_to?(:to_hash)

      @upload = Metis::Upload[:file_id=> @file.to_hash[:id]]
    end
  end

  # Check the first kilobyte hash.
  def first_kilobyte_check

    raise_err(:BAD_REQ, 0, __method__) if !@params.key?('first_kib_hash')
    first_kib_hash = nil
    File.open(derive_directory+'/'+@params['file_name']+'.part') do |file|

      first_kib_hash = Digest::MD5.hexdigest(file.read(1024))
    end
    raise_err(:BAD_REQ, 0, __method__) if !first_kib_hash

    if first_kib_hash != @params['first_kib_hash']

      raise_err(:BAD_REQ, 0, __method__)
    end
  end

  # Extra items from the server will be added to these when the HMAC is
  # generated.
  def has_auth_params?(params)

    has_params = true
    auth_params = [

      'original_name',
      'file_name',
      'file_size',
      'group_name',
      'project_name',
      'project_name_full'
    ]

    auth_params.each do |auth_param|

      if !params.key?(auth_param) then has_params = false end
    end
    return has_params
  end

  # Add the extra items needed to generate an HMAC auth.
  def add_extra_auth_params

    @params['hashing_algorithm'] = 'MD5'
    @params['start_timestamp'] = Time::now.to_i
    @params['token'] = @user.token
    @params['user_email'] = @user.email
  end

  # Generate the HMAC.
  def generate_hmac
    ordered_params = SignService::order_params(@params)
    hmac_signature = SignService::sign_request(ordered_params, 'sha256')
  end

  def normalize_name(name)
    nm = name.dup
    nm.gsub!('/ ', '_')
    nm.gsub!('(', '')
    nm.gsub!(')', '')
    nm.gsub!(' ', '_')
    nm.gsub!('/', '_')
    nm.gsub!('-', '_')
    return nm
  end

  def directory_exists?
    Metis::File.directory?(derive_directory)
  end

  def file_exists?
    file_name = @params['file_name'].to_s
    Metis::File.file?(derive_directory+'/'+file_name)
  end

  def partial_exists?

    file_name = @params['file_name'].to_s
    File.file?(derive_directory+'/'+file_name+'.part')
  end

  def db_metadata_exists?

    file = FileModel::File[

      :group_name=> normalize_name(@params['group_name'].to_s),
      :project_name=> normalize_name(@params['project_name'].to_s),
      :file_name=> @params['file_name'].to_s
    ]
    return if file ? true : false
  end

  def derive_directory
    Secrets::ROOT_DIR+'/'+normalize_name(@params['project_name'].to_s)
  end

  def upload_params_valid?

    valid = true
    Conf::UPLOAD_VALIDATION_ITEMS.each do |key, value|

      if !@params.key?(key) then valid = false end
    end
    return valid
  end

  def file_params_valid?

    valid = true
    Conf::FILE_VALIDATION_ITEMS.each do |key, value|

      if !@params.key?(key) then valid = false end
    end
    return valid
  end

  # We need to add a check that if a key/value is not of the correct type, then
  # we need to remove it. Then when we use our 'upload_params_valid' and 
  # 'file_params_valid' the result will be false thus protecting us from invalid
  # data types.

  def update_upload_metadata
    params = @request.POST
    temp_file_path = @request['blob'][:tempfile].path
    partial_file_name = derive_directory+'/'+@params['file_name']+'.part'

    @upload.update(
      current_byte_position: File.size(partial_file_name),
      current_blob_size: File.size(temp_file_path),
      next_blob_size: @params['next_blob_size'],
      next_blob_hash: @params['next_blob_hash'],
      status: 'active'
    )
  end

  def remove_temp_file
    partial_file_name = derive_directory+'/'+@params['file_name']+'.part'
    if File.file?(partial_file_name) then File.delete(partial_file_name) end
  end

  def remove_file_on_disk
    file_name = derive_directory+'/'+@params['file_name']
    if File.file?(file_name) then File.delete(file_name) end
  end
end
