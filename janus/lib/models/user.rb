class User < Sequel::Model
  include Etna::Instrumentation

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

    payload = payload.symbolize_keys.except(:exp)

    user = User[email: payload[:email]]

    return nil unless user and payload == user.jwt_payload

    return user
  end

  def to_hash
    {
      email: email,
      name: name,
      flags: flags,
      joined_at: created_at.iso8601,
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

  def create_task_token!(project_name, **kwds)
    token_builder.create_task_token!(project_name, **kwds)
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

  time_it(:valid_signature?)

  def authorized?(pass)
    # A password can be 'nil' if one logs in via Shibboleth/MyAccess.
    return false unless pass_hash

    client_hash = Janus.instance.sign.hash_password(pass)
    return pass_hash == client_hash
  end

  def superuser?
    @superuser ||= has_role?('administration', 'administrator')
  end

  def supereditor?
    @supereditor ||= has_role?('administration', 'administrator', 'editor')
  end

  def superviewer?
    @superviewer ||= has_role?('administration', 'administrator', 'editor', 'viewer')
  end

  def has_role?(project,*roles)
    permissions.any? do |permission|
      roles.include?(permission.role) &&
        permission.project&.project_name == project
    end
  end

  def set_guest_permissions!(project_name)
    # For given project: if agreed, ensure at least guest access; if explicitly !agreed, remove guest access; if no cc aggreement, do nothing
    cc_agreements_project = CcAgreement.where(user_email: email, project_name: project_name).all
    most_recent = cc_agreements_project.sort_by(&:created_at).last

    project = Project[project_name: project_name]
    if project.nil?
      return failure(404, "Project #{project_name} does not exist.")
    end

    perm = Permission.where(user: self, project: project).first
    if most_recent
      if most_recent.agreed
        if !perm
          Permission.create(project: project, user: self, role: 'guest')
        end
      end
      if !most_recent.agreed
        if perm && perm.role == 'guest'
          perm.delete
        end
      end
    end

    self.permissions(reload: true)
  end
end
