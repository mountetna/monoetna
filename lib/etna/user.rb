module Etna
  class User
    attr_reader :first, :last, :email
    def initialize params
      @first, @last, @email, perm = params.values_at(:first, :last, :email, :perm)
      raise ArgumentError, "No email given!" unless @email
      @roles = Hash[
        perm.split(/\;/).map do |e|
          role, projects = e.split(/:/)
          projects.split(/\,/).map do |p|
            [ p, role ]
          end
        end.inject([], &:+)
      ]
    end

    ROLE_MATCH = {
      admin: /[Aa]/,
      editor: /[Ee]/,
      viewer: /[Vv]/,
      restricted: /[AEV]/,
    }
    def has_roles(project, *roles)
      role = @roles[project.to_s]

      return false unless role

      return roles.any? { |role_name| role.match(ROLE_MATCH[role_name]) }
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
      has_roles(project, :restricted)
    end

    def is_admin? project
      is_superuser? || has_roles(project, :admin)
    end
  end
end
