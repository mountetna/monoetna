class Metis
  class File < Sequel::Model
    def create_partial!
      MetisFile.new(partial_file_name, 'w').close
    end

    def partial_file_name
      MetisFile.join(project_dir, file_name + '.part')
    end

    def project_dir
      MetisFile.join(Metis.instance.config(:root_dir), project_name)
    end
  end
end
