require_relative '../nonce'

class AuthorizationController < Janus::Controller
  def login
    # Make sure the refer url is valid.
    unless refer_valid?(@params[:refer])
      raise Etna::BadRequest, 'Invalid url refer'
    end

    token = @request.cookies[Janus.instance.config(:token_name)]

    begin
      user = User.from_token(token)

      raise unless user

      # they have a valid token
      @response.redirect(@params[:refer], 302)
      @response.finish
    rescue
      if Janus.instance.config(:auth_method) == 'shibboleth'
        return login_shib
      else
        return login_form
      end
    end
  end

  def validate_login
    require_params(:email, :password)

    unless email_valid?(@params[:email])
     raise Etna::BadRequest, 'Invalid email'
    end

    # Make sure the refer url is valid.
    unless refer_valid?(@params[:refer])
      raise Etna::BadRequest, 'Invalid url refer'
    end

    # Get and check user and then check the password.
    user = User[email: @params[:email]]

    unless user && user.authorized?(@params[:password])
      raise Etna::BadRequest, 'Invalid login'
    end

    # Create a new token for the user.
    token = user.create_token!

    # On success return the user info.
    respond_with_cookie(token, @params[:refer])
  end

  # generates a nonce for users to sign
  def time_signature
    # The message is the current moment in ISO
    # format
    return success(
      Janus::Nonce.new(DateTime.now.iso8601).to_s
    )
  end

  def generate
    # The user must authorize the request either
    # 1. with their signature

    signed_nonce = (@request.env['HTTP_AUTHORIZATION'] || '')[ /\ASigned-Nonce (.*)\z/, 1 ]

    # 2. with their token
    token = (@request.env['HTTP_AUTHORIZATION'] || '')[ /\AEtna (.*)\z/, 1 ]

    token_type = @params[:token_type] || 'login'

    user = @user ? User[email: @user.email] : signed_nonce ? User.from_signed_nonce(signed_nonce) : User.from_token(token)

    raise Etna::Unauthorized, user if user.is_a?(String)

    if token_type == 'task'
      raise Etna::BadRequest, "No project_name specified!" unless @params[:project_name]

      begin
        return success(user.create_task_token!(@params[:project_name], read_only: @params[:read_only]))
      rescue Token::Error => e
        Janus.instance.logger.log_error(e)
        raise Etna::Unauthorized, "failed to create token"
      end
    end

    return success(user.create_token!)
  end

  def build
    require_params(:email, :perm, :exp)

    @janus_user = User[email: @params[:email] ]

    raise Etna::BadRequest, 'User not found' unless @janus_user

    token = Janus.instance.sign.jwt_token(
      {
        email: @janus_user.email,
        name: @janus_user.name,
        perm: @params[:perm],
        flags: @janus_user.flags&.join(';'),
        exp: @params[:exp].to_i
      }
    )
    
    checker = Token::Checker.new(token)

    # forbid bad permissions
    raise Etna::BadRequest, "Cannot make a token with invalid permissions!" unless checker.valid_permissions?

    raise Etna::BadRequest, "Token is invalid!" unless checker.valid_token?

    success(token)
  end

  def validate_task
    @janus_user = User[email: @user.email]

    raise Etna::Forbidden, 'User not found' unless @janus_user

    raise Etna::Unauthorized, 'Invalid task token' unless Token::Checker.new(@user.token).valid_task_token?(@janus_user)

    success_json(success: true)
  end

  private

  def noauth
    failure(401, 'You are unauthorized')
  end

  def login_shib
    # Check that this request came from shibboleth(shibd).
    email = (@request.env['HTTP_X_SHIB_ATTRIBUTE'] || '').downcase
    raise Etna::Unauthorized if email == '(null)' || email.empty?

    # Get and check user. No password required.
    user = User[email: email]
    raise Etna::Unauthorized, "Unauthorized request for #{email}" unless user

    # Create a new token for the user.
    token = user.create_token!

    respond_with_cookie(token, @params[:refer])
  end

  def login_form
    erb_view(:login_form)
  end

  def respond_with_cookie(token, refer)
    # Set redirect.
    @response.redirect(refer.gsub("http://", "https://"), 302)

    Janus.instance.set_token_cookie(@response,token)
    return @response.finish
  end

  private

  def refer_valid?(refer)
    # Attempt to parse the refer.
    uri = URI.parse(refer)
    host = uri.host.split('.')[-2,2].join('.') # Extract the root host.

    # Check to make sure the refer comes from the same domain as the token.
    host == Janus.instance.config(:token_domain) ? true : false
  rescue
    return false
  end

  # Quick check that the email is in a somewhat valid format.
  def email_valid?(eml)
    eml =~ /\A[\w+\-.]+@[a-z\d\-]+(\.[a-z\d\-]+)*\.[a-z]+\z/
  end

end
