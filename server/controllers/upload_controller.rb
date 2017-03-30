class UploadController < BasicController

  def run()

    set_user()
    raise_err(:BAD_REQ, 2, __method__) if !@user.valid?()

    return send(@action).to_json()
  end

  def authorize_upload()

    m = __method__

    # Check that the correct parameters are present.
    raise_err(:BAD_REQ, 0, m) if !has_auth_params?(@params)

    group_name = normalize_name(@params['group_name'].to_s())
    project_name = normalize_name(@params['project_name'].to_s())
    directory = Conf::ROOT_DIR+'/'+group_name+'/'+project_name

    # Check that the user has permission on project requested.
    if !@user.project_editor?(@params['project_name']) &&
       !@user.project_admin?(@params['project_name'])

      raise_err(:BAD_REQ, 3, m)
    end

    # Check that this file system is in sync with the auth server. If there is a
    # group/project set in Janus there should be a corresponding directory
    # in Metis.
    raise_err(:SERVER_ERR,0,m) if !File.directory?(directory)

    # Generate a File Model and check to see if it already has metadata in the
    # system or a file on the disk.
    @file = FileModel.new(@params, @redis_service)
    if @file.file_exists?() || @file.db_metadata_exists?()

      raise_err(:BAD_REQ,5,m)
    end

    return generate_hmac_authorization()
  end

  # start_upload will create a metadata entry in the database and also a file on
  # the file system with 0 bytes.
  def start_upload()

    m = __method__

    # Check that the initialize/start request is valid.
    if !hmac_valid?() then return end

    @file = FileModel.new(@params, @redis_service)

    # Check that the upload directory exists. It should.
    raise_err(:SERVER_ERR, 0, m) if !@file.directory_exists?()

    # Make sure that a previous entry does not exist in the database or on the 
    # filesystem.
    raise_err(:BAD_REQ, 5, m) if @file.file_exists?()
    raise_err(:BAD_REQ, 6, m) if @file.db_metadata_exists?()

    @file.set_file_on_system!()

    # Make sure that everything was set ok.
    #raise_err(:BAD_REQ, 10, m) if !@file.file_exists?()
    #raise_err(:SERVER_ERR, 1, m) if !@file.db_metadata_exists?()

    response = {

      :success=> true,
      :request=> nil,
      :signature=> nil,
      :byte_count=> nil,
      :status=> 'initialized'
    }
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
      return {:success=> false, :request=> @params, :status=> 'not authorized'}
    end

    @params['status'] = 'authorized'
    return {

      :success=> true,
      :request=> @params,
      :hmac_signature=> generate_hmac(),
      :status=> 'authorized'
    }
  end

  # Add the extra items needed to generate an HMAC auth.
  def add_extra_auth_params()

    @params['directory'] = @file.directory
    @params['hashing_algorithm'] = 'MD5'
    @params['start_timestamp'] = Time::now.to_i
    @params['token'] = @user.token
    @params['user_email'] = @user.email
    @params['user_id'] = @user.id
    @params['old_index'] = @params['db_index']
    @params['db_index'] = @redis_service.get_new_index()
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
end