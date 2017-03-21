class UploadController < BasicController

  def run()

    set_user()
    if !@user.valid?() then raise_err(:BAD_REQ, 2, __method__) end

    return send(@action).to_json()
  end

  def authorize_upload()

    m = __method__

    # Check that the correct parameters are present.
    if !has_auth_params?(@params) then raise_err(:BAD_REQ, 0, m) end

    group_id = @params['group_id']
    project_id = @params['project_id']
    directory = Conf::ROOT_DIR+'/'+group_id.to_s+'/'+project_id.to_s()

    # Check that the user has permission on project requested.
    if !@user.project_editor?(project_id) && !@user.project_admin?(project_id)

      raise_err(:BAD_REQ, 3, m)
    end

    # Check that this file system is in sync with the auth server. If there is a
    # group/project set in Janus there should be a corresponding directory
    # in Metis.
    if !File.directory?(directory) then raise_err(:SERVER_ERR,0,m) end

    # Generate a File Model and check to see if it already has metadata in the
    # system or a file on the disk.
    @file = FileModel.new(@params, @redis_service)
    if @file.file_exists?() || @file.metadata_exists?()

      raise_err(:BAD_REQ, 5, m)
    end

    return generate_hmac_authorization()
  end

  private
  def has_auth_params?(params)

    has_params = true
    auth_params = [

      'token',
      'project_name',
      'project_id',
      'role',
      'file_name',
      'db_index',
      'group_id'
    ]

    auth_params.each do |auth_param|

      if !params.key?(auth_param) then has_params = false end
    end
    return has_params
  end

  def generate_hmac_authorization()

    # Add the extra items needed to generate an HMAC auth.
    @params['directory'] = @file.directory
    @params['expires'] = Conf::UPLOAD_EXPIRE
    @params['signing_algorithm'] = 'sha256'
    @params['hashing_algorithm'] = 'MD5'
    @params['start_timestamp'] = Time::now.to_i
    @params['token'] = @user.token
    @params['user_email'] = @user.email
    @params['user_id'] = @user.id
    @params['old_index'] = @params['db_index']
    @params['db_index'] = @redis_service.get_new_index()

    # Generate the HMAC.
    ordered_params = SignService::order_params(@params)
    hmac_signature = SignService::sign_request(ordered_params, 'sha256')
    @params['status'] = 'authorized'

    # Prep and send the response.
    response = {

      :success=> true,
      :request=> @params,
      :hmac_signature=> hmac_signature,
      :status=> 'authorized'
    }
  end
end