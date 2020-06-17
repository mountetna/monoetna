class UserController < Janus::Controller
  def update_key
    @janus_user = User[email: @user.email]

    raise Etna::Forbidden, 'User not found' unless @janus_user

    # parse the key
    begin
      key = OpenSSL::PKey.read(@params[:pem])
      raise Etna::BadRequest, 'Key must be RSA' unless key.is_a?(OpenSSL::PKey::RSA)
      raise Etna::BadRequest, 'Cannot use a private key' if key.private?
      raise Etna::BadRequest, 'Key must be 2048 bits' unless key.n.num_bits >= 2048
    rescue OpenSSL::PKey::PKeyError
      raise Etna::BadRequest, 'Improperly formatted key'
    end

    @janus_user.public_key = @params[:pem]
    @janus_user.save

    success('User Public Key updated')
  end
end
