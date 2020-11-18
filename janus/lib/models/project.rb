class Project < Sequel::Model
  one_to_many :permissions
  many_to_one :group

  PROJECT_NAME_MATCH=/\A[a-z][a-z0-9]*\Z/

  def self.valid_name?(project_name)
    project_name =~ PROJECT_NAME_MATCH && !project_name.start_with?('pg_')
  end

  def to_hash
    {
      project_id: id,
      group_id: group_id,
      group_name: group.group_name,
      project_name: project_name,
      project_name_full: project_name_full
    }
  end
end
