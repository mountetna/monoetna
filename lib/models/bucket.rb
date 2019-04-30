class Metis
  class Bucket < Sequel::Model
    one_to_many :files

    def allowed?(user)
      name == 'files'
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
