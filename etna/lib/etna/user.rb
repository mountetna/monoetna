require_relative './permissions'

module Etna
  class User
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
      @permissions ||= Etna::Permissions.from_encoded_permissions(@encoded_permissions).to_hash
    end

    def has_flag?(flag)
      @flags.include?(flag)
    end

    def projects
      permissions.keys
    end

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
      is_superviewer? || has_any_role?(project, :admin, :editor, :viewer)
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
  end
end
