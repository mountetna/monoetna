class User < Sequel::Model
  one_to_many :permissions

  def validate
    super
    errors.add(:email, 'must be lowercase') if email =~ /[A-Z]/
  end

  def self.from_signed_nonce(signed_nonce)
    nonce, email, signature = signed_nonce.split(/\./).map.with_index do |p,i|
      i == 0 ? p : Base64.decode64(p)
    end

    # validate the nonce
    return "invalid nonce #{nonce}" unless Janus::Nonce.valid_nonce?(nonce)

    # validate the email
    return "invalid email" if email =~ /[^[:print:]]/

    # find the user
    user = User[email: email]

    return "no user" unless user

    # check the user's signature
    txt_to_sign = signed_nonce.split('.')[0..1].join('.')

    return "invalid signature" unless user.valid_signature?(txt_to_sign, signature)

    return user
  end

  def self.from_token(token)
    payload, header = Janus.instance.sign.jwt_decode(token)

    payload = payload.transform_keys(&:to_sym)
    payload.delete(:exp)

    user = User[email: payload[:email]]

    return nil unless user and payload == user.jwt_payload

    return user
  end

  def to_hash
    {
      email: email,
      name: name,
      flags: flags,
      public_key: public_key && key_fingerprint
    }.compact
  end

  def key_fingerprint
    pkey = OpenSSL::PKey::RSA.new(public_key)

    data_string = [7].pack('N') + 'ssh-rsa' + pkey.public_key.e.to_s(0) + pkey.public_key.n.to_s(0)

    OpenSSL::Digest::MD5.hexdigest(data_string).scan(/../).join(':')
  end

  def token_builder
    @token_builder ||= Token::Builder.new(self)
  end

  def jwt_payload
    token_builder.jwt_payload
  end

  def create_token!
    token_builder.create_token!
  end

  def create_task_token!(project_name)
    token_builder.create_task_token!(project_name)
  end

  def create_viewer_token!
    token_builder.create_viewer_token!
  end

  def valid_signature?(text, signature)
    return nil unless public_key

    pkey = OpenSSL::PKey::RSA.new(public_key)

    verified = pkey.verify(
      OpenSSL::Digest::SHA256.new,
      signature, text
    )
    OpenSSL.errors.clear

    return verified
  end

  def authorized?(pass)
    # A password can be 'nil' if one logs in via Shibboleth/MyAccess.
    return false unless pass_hash

    client_hash = Janus.instance.sign.hash_password(pass)
    return pass_hash == client_hash
  end

  def superuser?
    @superuser ||= permissions.any? do |permission|
      permission.role == 'administrator' &&
        permission.project &&
        permission.project.project_name == 'administration'
    end
  end
end
