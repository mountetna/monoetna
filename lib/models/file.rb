class Metis
  class File < Sequel::Model
    def create_partial!
      ::File.new(partial_file_name, 'w').close
    end

    def partial_file_name
      ::File.join(project_dir, file_name + '.part')
    end
    
    def has_data?
      file_name && ::File.exists?(location)
    end

    def location
      ::File.expand_path(::File.join(Metis.instance.project_dir(project_name), file_name))
    end
  end
end
