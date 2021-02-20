class Project < Sequel::Model
  one_to_many :permissions

  PROJECT_NAME_MATCH=/\A[a-z][a-z0-9]*(_[a-z][a-z0-9]*)*\Z/

  def self.valid_name?(project_name)
    project_name =~ PROJECT_NAME_MATCH && !project_name.start_with?('pg_')
  end
end
