class Metis
  class Bucket < Sequel::Model
    one_to_many :files

    ROLES = {
      # stupidly, Etna::User calls this 'admin'
      admin: 1,
      administrator: 1,
      editor: 2,
      viewer: 3
    }

    def allowed?(user)
      # Admins can always see the bucket
      return true if user.is_admin?(project_name)

      # Access is restricted by role
      if ROLES[access.to_sym]
        return ROLES[access.to_sym] >= ROLES[user.permissions[project_name][:role]]
      end

      # Access is restricted to a user list by email id
      if access =~ /@/
        return access.split(/,/).include?(user.email)
      end

      return false
    end

    def location
     ::File.expand_path(::File.join(
       Metis.instance.project_path(project_name),
       name
     ))
    end

    def to_hash
      {
        bucket_name: name,
        project_name: project_name,
        count: files.count
      }
    end
  end
end
