class Metis
  class File < Sequel::Model
    def create_partial!
      Metis::File.new(partial_file_name, 'w').close
    end

    def partial_file_name
      Metis::File.join(project_dir, file_name + '.part')
    end

    def project_dir
      Metis::File.join(Metis.instance.config(:root_dir), project_name)
    end
  end
end
