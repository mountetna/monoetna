class AdminController < Janus::Controller
  def main
    @janus_user = User[email: @user.email]
    erb_view(:main)
  end

  def project
    @project = Project[project_name: @params[:project_name]]
    @static = nil
    if @user.is_superuser?
      @roles = [ 'administrator', 'viewer', 'editor', 'disabled' ]
    elsif @user.is_admin?(@params[:project_name])
      @roles = [ 'viewer', 'editor', 'disabled' ]
    else
      @roles = []
      @static = true
    end

    @project_roles = @project.permissions.group_by(&:role)
    erb_view(:project)
  end

  def update_permission
    require_param(:email)
    @project = Project[project_name: @params[:project_name]]

    permission = @project.permissions.find do |p| p.user.email == @params[:email] end

    # fix strings from HTML POST
    @params[:privileged] = @params[:privileged] == 'true' if [ 'true', 'false' ].include?(@params[:privileged])

    raise Etna::BadRequest, "No such user on project #{@params[:project_name]}!" unless permission

    raise Etna::Forbidden, 'Cannot update admin role!' if permission.role == 'administrator' && !@user.is_superuser?

    raise Etna::Forbidden, 'Cannot grant admin role!' if @params[:role] == 'administrator' && !@user.is_superuser?

    raise Etna::BadRequest, "Unknown role #{@params[:role]}" unless !@params[:role] || [ 'administrator', 'viewer', 'editor', 'disabled' ].include?(@params[:role])

    if @params[:role] == 'disabled'
      permission.delete
    else
      permission.role = @params[:role] if @params[:role]
      permission.privileged = @params[:privileged] if [true,false].include?(@params[:privileged])
      permission.affiliation = @params[:affiliation] if @params[:affiliation]
      permission.save
    end

    @response.redirect("/project/#{@params[:project_name]}")
    @response.finish
  end

  def add_user
    require_params(:email, :name, :role)
    @project = Project[project_name: @params[:project_name]]

    raise Etna::Forbidden, 'Cannot set admin role!' if @params[:role] == 'administrator'

    raise Etna::BadRequest, "Unknown role #{@params[:role]}" unless [ 'viewer', 'editor' ].include?(@params[:role])

    if @project.permissions.any? { |p| p.user.email == @params[:email] }
      raise Etna::BadRequest, "Duplicate permission on project #{@params[:project_name]}!"
    end

    user = User[email: @params[:email]]
    unless user
      raise Etna::BadRequest, 'Badly formed email address' unless @params[:email] =~ URI::MailTo::EMAIL_REGEXP

      names = @params[:name].split
      raise Etna::BadRequest, 'Missing name' if names.empty?
      first, last = names.length > 1 ? [ names[0..-2].join(' '), names.last ] : names

      user = User.create(email: @params[:email], first_name: first, last_name: last)
    end

    permission = Permission.create(project: @project, user: user, role: @params[:role])
    permission.role = @params[:role] if [ 'viewer', 'editor' ].include?(@params[:role])
    permission.privileged = false
    permission.affiliation = @params[:affiliation]
    permission.save

    @response.redirect("/project/#{@params[:project_name]}")
    @response.finish
  end

  def add_project
    require_params(:project_name, :project_name_full)

    raise Etna::BadRequest, "project_name should be like #{Project::PROJECT_NAME_MATCH.source}" unless Project.valid_name?(@params[:project_name])

    raise Etna::BadRequest, 'project_name_full cannot be empty' if @params[:project_name_full].nil? || @params[:project_name_full].empty?

    project = Project[project_name: @params[:project_name]]

    raise Etna::BadRequest, 'Duplicate project_name' if project

    project = Project.create(
      project_name: @params[:project_name],
      project_name_full: @params[:project_name_full]
    )

    @response.redirect('/')
    @response.finish
  end
end
