class UploadController < BasicController

  def run()

    set_user()
    raise_err(:BAD_REQ, 3, __method__) if !@user.valid?()
    normalize_params()
    return send(@action).to_json()
  end

  def authorize_upload()

    auth_upload_error_check()
    return generate_hmac_authorization()
  end

  # start_upload will create a metadata entry in the database and also a file on
  # the file system with 0 bytes.
  def start_upload()

    start_upload_error_check()

    # Change the status of the upload.
    @params['status'] = 'initialized'

    # Write the metadata to postgres.
    @file = PostgresService.create_new_file!(@params)
    @upload = PostgresService.create_new_upload!(@params, @file.to_hash[:id])

    # Create a partial file on the system.
    create_partial_file()

    # Send upload initiated
    return { :success=> true, :request=> @params }
  end

  # Upload a chunk of the file.
  def upload_blob()

    upload_blob_error_check()

    append_blob()
    update_upload_metadata()

    if upload_complete?()

      make_file_permanent()
      return send_upload_complete()
    else

      return send_upload_status()
    end
  end

  def pause_upload()

    common_error_check()
    @upload.update(:status=> 'paused')
    return send_upload_status()
  end

  def cancel_upload()

    common_error_check()
    @upload.delete() # Remove metadata
    @file.delete()   # Remove metadata
    remove_temp_file()
    @params.delete('blob')
    @params['status'] = 'cancelled'
    return { :success=> true, :request=> @params }
  end

  def remove_file()

    remove_file_error_check()
    @file.delete()   # Remove metadata
    remove_file_on_disk()
    @params['status'] = 'removed'
    return { :success=> true, :request=> @params }
  end

  def remove_failed()

    remove_failed_error_check()
    @upload.delete() # Remove metadata
    @file.delete()   # Remove metadata
    remove_temp_file()
    @params['status'] = 'removed'
    return { :success=> true, :request=> @params }
  end

  def recover_upload()

    recover_upload_error_check()

    raise_err(:BAD_REQ, 0, __method__) if !@params.key?('last_kib_hash')

    @upload.update(

      :status=> 'paused',
      :next_blob_size=> 1024,
      :next_blob_hash=> @params['last_kib_hash']
    )

    @params['hashing_algorithm'] = 'MD5'
    @params['start_timestamp'] = @file.to_hash()[:start_upload].to_time().to_i()
    @params['token'] = @user.token
    @params['user_email'] = @user.email
    @params['hmac_signature'] = generate_hmac()

    @params.delete('first_kib_hash')
    @params.delete('last_kib_hash')
    @params['status'] = 'paused'
    return { :success=> true, :request=> @params }
  end

  private

  # All the items to check before we authorize an upload.
  def auth_upload_error_check()

    # Check that the correct parameters are present.
    raise_err(:BAD_REQ, 0, __method__) if !has_auth_params?(@params)

    # Check that the user has edit permission on the project requested.
    user_edit_check()

    # Check that this file system is in sync with the auth server. If there is a
    # group/project set in Janus there should be a corresponding directory
    # in Metis.
    raise_err(:SERVER_ERR, 2, __method__) if !directory_exists?()

    # Check that the file does not exist on disk or in the db.
    raise_err(:BAD_REQ, 4, __method__) if file_exists?() || partial_exists?()
    raise_err(:BAD_REQ, 5, __method__) if db_metadata_exists?()
  end

  # Details to check before we start an upload.
  def start_upload_error_check()

    # Check the HMAC
    raise_err(:BAD_REQ, 8, __method__) if !hmac_valid?()

    # Check that the upload directory exists. It should.
    raise_err(:SERVER_ERR, 2, __method__) if !directory_exists?()

    # Check that the file does not exist on disk or in the db.
    raise_err(:BAD_REQ, 4, __method__) if file_exists?() || partial_exists?()
    raise_err(:BAD_REQ, 5, __method__) if db_metadata_exists?()

    # Check the params from the client.
    if !upload_params_valid?() || !file_params_valid?()

      raise_err(:BAD_REQ, 0, __method__)
    end
  end

  # Things to check before we upload a blob.
  def upload_blob_error_check()

    common_error_check()

    # Check that the partial file exists.
    raise_err(:BAD_REQ, 6, __method__) if !partial_exists?()

    # Check that the full file DOES NOT exist.
    raise_err(:BAD_REQ, 4, __method__) if file_exists?()

    # Check the uploaded blob's integrity.
    next_blob_size = @upload.to_hash[:next_blob_size]
    next_blob_hash = @upload.to_hash[:next_blob_hash]

    if !blob_integrity_ok?(next_blob_size, next_blob_hash)

      raise_err(:BAD_REQ, 9, __method__) 
    end
  end

  # Some common things we should check before we make an action.
  def common_error_check()

    # Check the HMAC
    raise_err(:BAD_REQ, 8, __method__) if !hmac_valid?()

    # Check that the file metadata exists and set global vars if they exist.
    file_metadata_check_and_set()
  end

  # All the things we should check before we remove a file.
  def remove_file_error_check()

    # Check that the user has edit permission on the project requested.
    user_edit_check()

    # Check that the file metadata exists.
    @file = FileModel::File[

      :group_name=> @params['group_name'],
      :project_name=> @params['project_name'],
      :file_name=> @params['file_name']
    ]
    raise_err(:SERVER_ERR, 7, __method__) if !@file

    # The 'upload' portion of the metadata should NOT exist.
    @upload = FileModel::Upload[:file_id=> @file.to_hash[:id]]
    raise_err(:SERVER_ERR, 5, __method__) if @upload

    # The file should exist on disk BUT the partial should NOT.
    raise_err(:BAD_REQ, 6, __method__) if !file_exists?()
    raise_err(:BAD_REQ, 4, __method__) if partial_exists?()
  end

  def remove_failed_error_check()

    # Check that the file metadata exists and set global vars if they exist.
    file_metadata_check_and_set()

    # The file should NOT exist but the partial should.
    raise_err(:BAD_REQ, 4, __method__) if file_exists?()
    raise_err(:BAD_REQ, 6, __method__) if !partial_exists?()

    # Check that the user requesting the 'remove' is the same as the uploader.
    raise_err(:BAD_REQ, 8, __method__) if !@user.email() == @file[:upload_by]
  end

  # All the items to check before we 'recover' an upload sequence.
  def recover_upload_error_check()

    # Check that the correct parameters are present.
    raise_err(:BAD_REQ, 0, __method__) if !has_auth_params?(@params)

    # Check that the user has edit permission on the project requested.
    user_edit_check()

    # Check that this file system is in sync with the auth server. If there is a
    # group/project set in Janus there should be a corresponding directory
    # in Metis.
    raise_err(:SERVER_ERR, 2, __method__) if !directory_exists?()

    # Check that a full file does not exist. A partial should exist.
    raise_err(:BAD_REQ, 4, __method__) if file_exists?()
    raise_err(:BAD_REQ, 6, __method__) if !partial_exists?()

    # Check that the file metadata exists and set global vars if they exist.
    file_metadata_check_and_set()

    # Check the first kilobyte hash.
    first_kilobyte_check()
  end

  # Check that the file metadata exists and set global vars if they exist.
  # This method also sets the variables @file and @update. We need to find a
  # cleaner way to initialize these variables.
  def file_metadata_check_and_set()

    @file = FileModel::File[

      :group_name=> @params['group_name'],
      :project_name=> @params['project_name'],
      :file_name=> @params['file_name']
    ]
    raise_err(:BAD_REQ, 7, __method__) if !@file

    # The 'upload' portion of the metadata should exist.
    @upload = FileModel::Upload[:file_id=> @file.to_hash[:id]]
    raise_err(:BAD_REQ, 7, __method__) if !@upload
  end

  # Check that the user has edit permission on the project requested.
  def user_edit_check()

    if !@user.project_editor?(@params['project_name']) &&
       !@user.project_admin?(@params['project_name'])

      raise_err(:BAD_REQ, 8, __method__)
    end
  end

  # Check the first kilobyte hash.
  def first_kilobyte_check()

    raise_err(:BAD_REQ, 0, __method__) if !@params.key?('first_kib_hash')
    first_kib_hash = nil
    File.open(derive_directory()+'/'+@params['file_name']+'.part') do |file|

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
      'project_name'
    ]

    auth_params.each do |auth_param|

      if !params.key?(auth_param) then has_params = false end
    end
    return has_params
  end

  def generate_hmac_authorization()

    add_extra_auth_params()

    if !hmac_params_valid?()

      @params['status'] = 'not authorized'
      return { :success=> false, :request=> @params }
    end

    @params['status'] = 'authorized'
    @params['hmac_signature'] = generate_hmac()
    return { :success=> true, :request=> @params }
  end

  # Add the extra items needed to generate an HMAC auth.
  def add_extra_auth_params()

    @params['hashing_algorithm'] = 'MD5'
    @params['start_timestamp'] = Time::now.to_i
    @params['token'] = @user.token
    @params['user_email'] = @user.email
  end

  def hmac_params_valid?()

    params_valid = true
    Conf::SIGNATURE_ITEMS.each do |item|

      if !@params.key?(item) then params_valid = false end
    end
    return params_valid
  end

  def hmac_valid?()

    if !@params.key?('hmac_signature') then return false end
    if !hmac_params_valid?() then return false end
    if @params['hmac_signature'] != generate_hmac() then return false end
    return true
  end

  # Generate the HMAC.
  def generate_hmac()

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

  def directory_exists?()

    File.directory?(derive_directory())
  end

  def file_exists?()

    file_name = @params['file_name'].to_s()
    File.file?(derive_directory()+'/'+file_name)
  end

  def partial_exists?()

    file_name = @params['file_name'].to_s()
    File.file?(derive_directory()+'/'+file_name+'.part')
  end

  def db_metadata_exists?()

    file = FileModel::File[

      :group_name=> normalize_name(@params['group_name'].to_s()),
      :project_name=> normalize_name(@params['project_name'].to_s()),
      :file_name=> @params['file_name'].to_s()
    ]
    return if file ? true : false
  end

  def derive_directory()

    group_name = normalize_name(@params['group_name'].to_s())
    project_name = normalize_name(@params['project_name'].to_s())
    return Conf::ROOT_DIR+'/'+group_name+'/'+project_name
  end

  def upload_params_valid?()

    valid = true
    Conf::UPLOAD_VALIDATION_ITEMS.each do |key, value|

      if !@params.key?(key) then valid = false end
    end
    return valid
  end

  def file_params_valid?()

    valid = true
    Conf::FILE_VALIDATION_ITEMS.each do |key, value|

      if !@params.key?(key) then valid = false end
    end
    return valid
  end

  # We need to add a check that if a key/value is not of the correct type, then
  # we need to remove it. Then when we use our 'upload_params_valid' and 'file_params_valid'
  # the result will be false thus protecting us from invalid data types.
  def normalize_params()

    @params.each do |key, value|

      if Conf::FILE_VALIDATION_ITEMS.key?(key)

        if Conf::FILE_VALIDATION_ITEMS[key] == Integer 

          @params[key] = @params[key].to_i
        end
      end

      if Conf::UPLOAD_VALIDATION_ITEMS.key?(key)

        if Conf::UPLOAD_VALIDATION_ITEMS[key] == Integer 

          @params[key] = @params[key].to_i
        end
      end
    end
  end

  def create_partial_file()

    partial_file_name = derive_directory()+'/'+@params['file_name']+'.part'
    partial_file = File.new(partial_file_name, 'w')
    partial_file.close()
  end

  def blob_integrity_ok?(next_blob_size, next_blob_hash)

    # Check the blob hash.
    if !blob_hash_ok?(next_blob_hash) then return false end

    # Check the blob size.
    temp_file_path = @request['blob'][:tempfile].path()
    if File.size(temp_file_path).to_i != next_blob_size then return false end

    return true
  end

  def blob_hash_ok?(md5_from_status)

    temp_file_path = @request['blob'][:tempfile].path()
    md5_of_temp_file = Digest::MD5.hexdigest(File.read(temp_file_path))
    return (md5_from_status == md5_of_temp_file) ? true : false
  end

  def append_blob()

    temp_file_name = @request['blob'][:tempfile].path()
    partial_file_name = derive_directory()+'/'+@params['file_name']+'.part'
    partial_file = File.open(partial_file_name, 'ab')
    temp_file = File.open(temp_file_name, 'rb')
    partial_file.write(temp_file.read())

    partial_file.close()
    temp_file.close()
  end

  def update_upload_metadata()

    params = @request.POST()
    temp_file_path = @request['blob'][:tempfile].path()
    partial_file_name = derive_directory()+'/'+@params['file_name']+'.part'

    @upload.update(

      :current_byte_position=> File.size(partial_file_name),
      :current_blob_size=> File.size(temp_file_path),
      :next_blob_size=> @params['next_blob_size'],
      :next_blob_hash=> @params['next_blob_hash'],
      :status=> 'active'
    )
  end

  def upload_complete?()

    file_size = @file.to_hash[:file_size]
    partial_file_name = derive_directory()+'/'+@params['file_name']+'.part'
    return File.size(partial_file_name) == file_size ? true : false
  end

  def make_file_permanent()

    # Rename the partial file.
    full_path = derive_directory()+'/'+@file.to_hash[:file_name]
    File.rename(full_path+'.part', full_path)

    # Generate the file hash.
    file_hash = Digest::MD5.hexdigest(File.read(full_path))

    # Update the postgres record.
    @file.update(:finish_upload=> Time::now, :hash=> file_hash)

    # Remove the upload record.
    @upload.delete
  end

  def send_upload_status()

    @upload.to_hash.each(){ |key, value| @params[key] = value }
    @params.delete(:file_id)
    @params.delete('blob')
    { :success=> true, :request=> @params }
  end

  def send_upload_complete()

    @file.to_hash.each(){ |key, value| @params[key] = value }
    @params.delete(:file_id)
    @params['status'] = 'complete'
    { :success=> true, :request=> @params }
  end

  def remove_temp_file()

    partial_file_name = derive_directory()+'/'+@params['file_name']+'.part'
    if File.file?(partial_file_name) then File.delete(partial_file_name) end
  end

  def remove_file_on_disk()

    file_name = derive_directory()+'/'+@params['file_name']
    if File.file?(file_name) then File.delete(file_name) end
  end
end
