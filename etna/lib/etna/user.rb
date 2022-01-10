module Etna
  class User
    ROLE_NAMES = {
      'A' => :admin,
      'E' => :editor,
      'V' => :viewer
    }

    def initialize params, token=nil
      @name, @email, @encoded_permissions, encoded_flags, @task = params.values_at(:name, :email, :perm, :flags, :task)

      @flags = encoded_flags&.split(/;/) || []
      @token = token unless !token
      raise ArgumentError, "No email given!" unless @email
    end

    attr_reader :name, :email, :token

    def task?
      !!@task
    end

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
    def has_any_role?(project, *roles)
      perm = permissions[project.to_s]

      return false unless perm

      return roles.map(&:to_sym).include?(perm[:role])
    end

    def is_superuser? project=nil
      has_any_role?(:administration, :admin)
    end

    def is_supereditor? project=nil
      has_any_role?(:administration, :admin, :editor)
    end

    def is_superviewer? project=nil
      has_any_role?(:administration, :admin, :editor, :viewer)
    end

    def can_edit? project
      is_supereditor? || has_any_role?(project, :admin, :editor)
    end

    def can_view? project
      is_superviewer? || has_any_role?(project, :admin, :editor, :viewer) || resource_project?(project)
    end

    # superusers - administrators of the Administration group - cannot
    # automatically see restricted data, they should be granted
    # project-specific access.
    def can_see_restricted? project
      perm = permissions[project.to_s]
      perm && perm[:restricted]
    end

    def is_admin? project
      is_superuser? || has_any_role?(project, :admin)
    end

    def active? project=nil
      permissions.keys.length > 0
    end

    def user_projects
      return [] unless has_janus_config?

      janus_client(token).get_projects.json[:projects]
    end

    def resource_project?(project_name)
      return false unless has_janus_config?

      project = user_projects.select { |p| p.project_name == project_name }.first
      
      return false unless project

      project[:resource]
    end

    def has_janus_config?
      application.config(:janus) && application.config(:janus)[:host]
    end

    def janus_client(token)
      Etna::Clients::Janus.new(
        token: token,
        host: application.config(:janus)[:host]
      )
    end

    def application
      @application ||= Etna::Application.instance
    end
  end
end
