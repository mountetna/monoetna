module Etna
  class User
    ROLE_NAMES = {
      'A' => :admin,
      'E' => :editor,
      'V' => :viewer
    }

    def initialize params, token=nil
      @name, @email, @encoded_permissions, encoded_flags = params.values_at(:name, :email, :perm, :flags)

      @flags = encoded_flags&.split(/;/) || []
      @token = token unless !token
      raise ArgumentError, "No email given!" unless @email
    end

    attr_reader :name, :email, :token

    def permissions
      @permissions ||= @encoded_permissions.split(/\;/).map do |roles|
        role, projects = roles.split(/:/)

        projects.split(/\,/).reduce([]) do |perms,project_name|
          perms.push([
            project_name,
            {
              role: ROLE_NAMES[role.upcase],
              restricted: role == role.upcase
            }
          ])
        end
      end.inject([],:+).to_h
    end

    def has_flag?(flag)
      @flags.include?(flag)
    end

    def projects
      permissions.keys
    end

    ROLE_MATCH = {
      admin: /[Aa]/,
      editor: /[Ee]/,
      viewer: /[Vv]/,
      restricted: /[AEV]/,
    }
    def has_roles(project, *roles)
      perm = permissions[project.to_s]

      return false unless perm

      return roles.map(&:to_sym).include?(perm[:role])
    end

    def is_superuser? project=nil
      has_roles(:administration, :admin)
    end

    def can_edit? project
      is_superuser? || has_roles(project, :admin, :editor)
    end

    def can_view? project
      is_superuser? || has_roles(project, :admin, :editor, :viewer)
    end

    # superusers - administrators of the Administration group - cannot
    # automatically see restricted data, they should be granted
    # project-specific access.
    def can_see_restricted? project
      perm = permissions[project.to_s]
      perm && perm[:restricted]
    end

    def is_admin? project
      is_superuser? || has_roles(project, :admin)
    end

    def active? project=nil
      permissions.keys.length > 0
    end
  end
end
