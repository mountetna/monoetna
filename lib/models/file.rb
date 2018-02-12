class Metis
  class File < Sequel::Model
    one_to_one :upload
    def self.has_file?(project_name, file_name)
      file = self.where(project_name: project_name, file_name: file_name).first

      return file && file.has_data?
    end

    def has_data?
      file_name && ::File.exists?(location)
    end

    def location
      ::File.expand_path(::File.join(Metis.instance.project_dir(project_name), file_name))
    end
  end
end
