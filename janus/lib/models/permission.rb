class Permission < Sequel::Model
  many_to_one :project
  many_to_one :user

  def to_hash
    {
      id: id,
      user_id: user_id,
      project_id: project_id,
      role: role,
      affiliation: affiliation,
      project_name: project.project_name,
      user_email: user.email,
      group_id: project.group_id,
      group_name: project.group.group_name
    }
  end

  alias_method :privileged?, :privileged

  def role_key
    privileged ? role[0].upcase : role[0]
  end

  def admin?
    role == 'administrator'
  end

  def editor?
    role == 'editor' || admin?
  end

  def project_role
    [ role_key, project.project_name ]
  end
end
