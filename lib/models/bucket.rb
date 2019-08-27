class Metis
  class Bucket < Sequel::Model
    plugin :timestamps, update_on_create: true

    one_to_many :files
    one_to_many :folders

    ROLES = {
      # stupidly, Etna::User calls this 'admin'
      admin: 1,
      administrator: 1,
      editor: 2,
      viewer: 3
    }

    def self.valid_bucket_name?(bucket_name)
      !!(bucket_name =~ /\A\w+\z/)
    end

    def self.valid_access?(access)
      ROLES.keys.include?(access.to_sym) ||
        access.split(/,/).all?{ |e| e.strip =~ URI::MailTo::EMAIL_REGEXP }
    end

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
       'buckets',
       name
     ))
    end

    def rename!(new_bucket_name)
      old_location = location

      update(name: new_bucket_name)
      refresh

      new_location = location
      FileUtils.mv(old_location, new_location)
    end

    def remove!
      ::Dir.delete(location)
      delete
    end

    def can_remove?
      has_directory? && Dir.entries(location).size <= 2
    end

    def create_actual_bucket!
      FileUtils.mkdir_p(location)
    end

    def has_directory?
      ::File.exists?(location) && ::File.directory?(location)
    end

    def to_hash
      {
        bucket_name: name,
        project_name: project_name,
        access: access,
        description: description,
        count: files.count
      }
    end
  end
end
