class User < Sequel::Model
  one_to_many :permissions

  def name
    [ first_name, last_name ].compact.join(' ')
  end

  def jwt_payload
    {
      email: email,
      first: first_name,
      last: last_name,

      # Encode permissions as a string e.g. "a:p1,p2;e:p3;v:p4"
      perm:  permissions.map(&:project_role).group_by(&:first)
        .sort_by(&:first).map do |role_key, project_roles|
        [ role_key, project_roles.map(&:last).sort.join(',') ].join(':')
      end.join(';')
    }
  end

  def key_fingerprint
    pkey = OpenSSL::PKey::RSA.new(public_key)

    data_string = [7].pack('N') + 'ssh-rsa' + pkey.public_key.e.to_s(0) + pkey.public_key.n.to_s(0)

    OpenSSL::Digest::MD5.hexdigest(data_string).scan(/../).join(':')
  end

  # WARNING! In the event of a shibboleth login 'pass_hash' == nil!
  def create_token!
    # Time is in seconds, nil = no expiration
    expires = Time.now.utc + Janus.instance.config(:token_life)

    Janus.instance.sign.jwt_token(
      jwt_payload.merge(exp: expires.to_i)
    )
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
