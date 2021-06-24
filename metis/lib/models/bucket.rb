class Metis
  class Bucket < Sequel::Model
    plugin :timestamps, update_on_create: true

    one_to_many :files
    one_to_many :folders
    one_to_many :uploads

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

    def allowed?(user, hmac)
      # if there is no user, or the bucket owner is not metis, it must be
      # hmac-signed by the bucket owner
      return hmac&.valid? && hmac.id == owner.to_sym if !user || owner != 'metis'

      # Admins can always see the bucket
      return true if user.is_admin?(project_name)

      # User not allowed if project not in their permissions.
      return false unless user.permissions[project_name]

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

    def rename!(new_bucket_name)
      update(name: new_bucket_name)
      refresh
    end

    def remove!
      # Stop all uploads
      uploads.each(&:delete_with_partial!)

      # delete the record
      delete
    end

    def can_remove?
      files.count == 0 && folders.count == 0
    end

    def to_hash
      {
        bucket_name: name,
        project_name: project_name,
        access: access,
        description: description,
        count: Metis::File.where(bucket: self).count
      }
    end
  end
end
