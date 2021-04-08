module Token
  class Error < StandardError
  end

  def self.role_key(role, privileged)
    privileged ? role[0].upcase : role[0]
  end

  ROLE_KEYS = {
    a: 'administrator',
    e: 'editor',
    v: 'viewer'
  }
  def self.role_from_key(role_key)
    return [ ROLE_KEYS[role_key.downcase.to_sym], role_key == role_key.upcase ]
  end

  class Builder
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
      filtered_perms.map do |permission|
        [ Token.role_key(permission.role, permission.privileged?), permission.project_name ]
      end.group_by(&:first)
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
        raise Token::Error, "Cannot write invalid permission on task token!"
      end

      # degrade admin permissions
      payload[:perm] = payload[:perm].tr('Aa', 'Ee') if payload[:perm] =~ /^[Aa]/

      # set task flag
      payload[:task] = true

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

  class Checker
    def initialize(token)
      @token = token
    end

    def valid_task_token?(janus_user)
      # the task flag is required
      return false unless payload[:task]

      # there must be only one project
      return false unless permissions.size == 1
      
      token_permission = permissions.first

      return false unless token_permission[:projects].size == 1

      project_name = token_permission[:projects].first
      token_role, token_privileged = Token.role_from_key(token_permission[:role])

      # no superuser or admin task tokens
      return false if project_name == 'administration'
      return false if token_role == 'administrator'

      # the project must be valid
      janus_permission = janus_user.permissions.find do |permission|
        permission.project_name == project_name
      end

      # the janus permission cannot be less than the token role
      return false unless janus_permission && janus_permission.role <= token_role

      return true
    end

    private
    def payload
      @payload ||= JSON.parse(parts[1], symbolize_names: true)
    end

    def header
      @header ||= JSON.parse(parts[0], symbolize_names: true)
    end

    def permissions
      @permissions ||= payload[:perm].split(';').map do |r|
        role, projects = r.split(':')
        {
          role: role,
          projects: projects.split(',')
        }
      end
    end

    def parts
      parts ||= @token.split('.').map do |p|
        Base64.decode64(p)
      end
    end
  end
end
