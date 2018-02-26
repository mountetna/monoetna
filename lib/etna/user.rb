module Etna
  class User
    ROLE_NAMES = {
      'A' => :admin,
      'E' => :editor,
      'V' => :viewer
    }

    attr_reader :first, :last, :email, :permissions
    def initialize params
      @first, @last, @email, perm = params.values_at(:first, :last, :email, :perm)
      raise ArgumentError, "No email given!" unless @email
      @permissions = Hash[
        perm.split(/\;/).map do |e|
          role, projects = e.split(/:/)
          projects.split(/\,/).map do |p|
            [ p, { role: ROLE_NAMES[role.upcase], restricted: role == role.upcase } ]
          end
        end.inject([], &:+)
      ]
    end

    def projects
      @permissions.keys
    end

    ROLE_MATCH = {
      admin: /[Aa]/,
      editor: /[Ee]/,
      viewer: /[Vv]/,
      restricted: /[AEV]/,
    }
    def has_roles(project, *roles)
      perm = @permissions[project.to_s]

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
      perm = @permissions[project.to_s]
      perm && perm[:restricted]
    end

    def is_admin? project
      is_superuser? || has_roles(project, :admin)
    end
  end
end
