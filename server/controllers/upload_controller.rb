class UploadController < BasicController

  def run()

    set_user()
    raise_err(:BAD_REQ, 2, __method__) if !@user.valid?()
    normalize_params()
    return send(@action).to_json()
  end

  def authorize_upload()

    # Check that the correct parameters are present.
    raise_err(:BAD_REQ, 0, __method__) if !has_auth_params?(@params)

    # Check that the user has permission on project requested.
    if !@user.project_editor?(@params['project_name']) &&
       !@user.project_admin?(@params['project_name'])

      raise_err(:BAD_REQ, 3, __method__)
    end

    # Check that this file system is in sync with the auth server. If there is a
    # group/project set in Janus there should be a corresponding directory
    # in Metis.
    raise_err(:SERVER_ERR, 0 , __method__) if !directory_exists?()

    # Check that the file does not exist on disk or in the db.
    if file_exists?() || partial_exists?() || db_metadata_exists?() 

      raise_err(:BAD_REQ, 5 , __method__)
    end 

    return generate_hmac_authorization()
  end

  # start_upload will create a metadata entry in the database and also a file on
  # the file system with 0 bytes.
  def start_upload()

    # Check the HMAC
    raise_err(:BAD_REQ, 7, __method__) if !hmac_valid?()

    # Check that the upload directory exists. It should.
    raise_err(:SERVER_ERR, 0, __method__) if !directory_exists?()

    # Check that the file does not exist on disk or in the db.
    if file_exists?() || partial_exists?() || db_metadata_exists?() 

      raise_err(:BAD_REQ, 5 , __method__)
    end

    # Check the params from the client.
    if !upload_params_valid?() || !file_params_valid?()

      raise_err(:BAD_REQ, 1, __method__)
    end

    # Set the data to the system.
    file = PostgresService.create_new_file!(@params).to_hash
    upload = PostgresService.create_new_upload!(@params, file[:id])
    return send_upload_initiated()
  end

  private
  def has_auth_params?(params)

    has_params = true
    auth_params = [

      'original_name',
      'file_name',
      'file_size',
      'group_id',
      'project_id'
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

  def send_upload_initiated()

    @params['hmac_signature'] = generate_hmac()
    @params['status'] = 'initialized'
    { :success=> true, :request=> @params }
  end

  # Generate the key used to access the file's metadata in Redis.
  def generate_metadata_key()
  
    items = ['db_index', 'group_name', 'project_name', 'file_name']
    items.map! do |item|
  
      if item == 'group_name' || item == 'project_name'
  
        normalize_name(@params[item])
      else
  
        @params[item]
      end
    end
    return items.join('-')
  end

  def blob_hash_ok?()

    md5_from_status = @file['next_blob_hash']
    puts md5_from_status
    temp_file_path = @request['blob'][:tempfile].path()
    puts temp_file_path
    md5_of_temp_file = Digest::MD5.hexdigest(File.read(temp_file_path))
    puts md5_of_temp_file

    return (md5_from_status == md5_of_temp_file) ? true : false
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
end