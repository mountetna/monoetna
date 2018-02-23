class Metis
  class File < Sequel::Model
    def self.md5(path)
      # use md5sum to avoid reading blob
      %x{ md5sum #{path} }.split.first
    end

    one_to_many :uploads
    def self.has_file?(project_name, file_name)
      file = self.where(project_name: project_name, file_name: file_name).first

      return file && file.has_data?
    end

    def has_data?
      file_name && ::File.exists?(location)
    end

    def location
      ::File.expand_path(::File.join(
        Metis.instance.project_path(project_name), file_name
      ))
    end

    def set_file_data(file_path)
      # Rename the existing file.
      ::File.rename(
        file_path,
        location
      )

      # update the hash
      self.update(hash: Metis::File.md5(location))
    end
  end
end
