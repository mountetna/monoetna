class Metis
  class Upload < Sequel::Model
    many_to_one :file

    def self.upload_started?(project_name, file_name)
      !self.where(project_name: project_name, file_name: file_name).empty?
    end
  end
end
