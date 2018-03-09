class Metis
  class File < Sequel::Model
    def self.md5(path)
      # use md5sum to avoid reading blob
      %x{ md5sum '#{path}' }.split.first
    end

    one_to_many :uploads
    def self.has_file?(project_name, file_name)
      file = self.where(project_name: project_name, file_name: file_name).first

      return file && file.has_data?
    end

    def compute_hash!
      if has_data?
        self.file_hash = Metis::File.md5(location)
        save
      end
    end

    def has_data?
      file_name && ::File.exists?(location)
    end

    def location
      ::File.expand_path(::File.join(
        Metis.instance.project_path(project_name), file_name
      ))
    end

    def to_hash
      {
        file_name: file_name,
        project_name: project_name,
        original_name: original_name,
        size: actual_size,
        file_hash: file_hash
      }
    end

    def actual_size
      has_data? ? ::File.size(location) : nil
    end

    def to_json(options)
      to_hash.to_json(options)
    end

    def set_file_data(file_path)
      # Rename the existing file.
      ::File.rename(
        file_path,
        location
      )

      # update the hash
      self.update(file_hash: Metis::File.md5(location))
    end
  end
end
