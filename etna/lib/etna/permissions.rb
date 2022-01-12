module Etna
  class Permissions
    def initialize(permissions)
      @permissions = permissions
    end

    def self.from_encoded_permissions(encoded_permissions)
      perms = encoded_permissions.split(/\;/).map do |roles|
        role, projects = roles.split(/:/)

        projects.split(/\,/).reduce([]) do |perms, project_name|
          perms << Etna::Permission.new(role, project_name)
        end
      end.flatten

      Etna::Permissions.new(perms)
    end

    def self.from_hash(permissions_hash)
      perms = permissions_hash.map do |project_name, role_hash|
        Etna::Permission.new(
          Etna::Role.new(role_hash[:role], role_hash[:restricted]).key,
          project_name
        )
      end

      Etna::Permissions.new(perms)
    end

    def to_string
      @permissions_string ||= @permissions.group_by(&:role_key)
        .sort_by(&:first)
        .map do |role_key, permissions|
        [
          role_key,
          permissions.map(&:project_name).sort.join(","),
        ].join(":")
      end.join(";")
    end

    def to_hash
      @permissions_hash ||= @permissions.map do |permission|
        [permission.project_name, permission.to_hash]
      end.to_h
    end

    def add_permission(permission)
      return if current_project_names.include?(permission.project_name)

      @permissions << permission
    end

    private

    def current_project_names
      @permissions.map(&:project_name)
    end
  end

  class Permission
    attr_reader :role, :project_name, :role_key

    ROLE_NAMES = {
      "A" => :admin,
      "E" => :editor,
      "V" => :viewer,
    }

    def initialize(role, project_name)
      @role_key = role
      @role = Etna::Role.new(ROLE_NAMES[role.upcase], role == role.upcase)
      @project_name = project_name
    end

    def to_hash
      role.to_hash
    end
  end

  class Role
    attr_reader :role, :restricted
    def initialize(role, restricted)
      @role = role
      @restricted = restricted
    end

    def key
      role_key = role.to_s[0]
      restricted ? role_key.upcase : role_key
    end

    def to_hash
      {
        role: role,
        restricted: restricted,
      }
    end
  end
end
