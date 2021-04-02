class TokenBuilder
  class Error < StandardError
  end

  def initialize(janus_user)
    @janus_user = janus_user
  end

  def jwt_payload(filters:{})
    {
      email: @janus_user.email,
      name: @janus_user.name,
      perm:  serialize_permissions(filters: filters),
      flags: @janus_user.flags&.join(';')
    }.compact
  end

  def serialize_permissions(filters: {})
    # Encode permissions as a string e.g. "a:p1,p2;e:p3;v:p4"
    filtered_perms = @janus_user.permissions.select do |perm|
      filters.all? do |k,v|
        perm.send(k) == v
      end
    end
    filtered_perms
      .map(&:project_role)
      .group_by(&:first)
      .sort_by(&:first)
      .map do |role_key, project_roles|
        [
          role_key,
          project_roles.map(&:last).sort.join(',')
        ].join(':')
    end.join(';')
  end

  def create_token!
    # Time is in seconds, nil = no expiration
    expires = Time.now.utc + Janus.instance.config(:token_life)

    Janus.instance.sign.jwt_token(
      jwt_payload.merge(exp: expires.to_i)
    )
  end

  def self.valid_task_token?(token, janus_user)
    # take apart the token
    payload = token.split('.')[1]
  end

  def create_task_token!(project_name)
    # Time is in seconds, nil = no expiration
    expires = Time.now.utc + Janus.instance.config(:task_token_life)

    payload = jwt_payload(filters: { project_name: project_name }).merge(exp: expires.to_i)

    unless payload[:perm] =~ /^[AaEeVv]:#{Project::PROJECT_NAME_MATCH.source}$/
      raise TokenBuilder::Error, "Cannot write invalid permission on task token!"
    end

    # degrade admin permissions
    payload[:perm] = payload[:perm].tr('Aa', 'Ee') if payload[:perm] =~ /^[Aa]/

    payload[:task] = 1

    Janus.instance.sign.jwt_token(
      payload
    )
  end

  def create_viewer_token!
    # Time is in seconds, nil = no expiration
    expires = Time.now.utc + Janus.instance.config(:token_life)

    Janus.instance.sign.jwt_token(
      jwt_payload(filters: { role: 'viewer'}).merge(exp: expires.to_i)
    )
  end

end
