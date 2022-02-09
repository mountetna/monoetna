class AdminController < Janus::Controller
  def main
    @janus_user = User[email: @user.email]
    erb_view(:client)
  end

  def project
    @project = Project[project_name: @params[:project_name]]

    raise Etna::BadRequest, "No such project #{@params[:project_name]}" unless @project

    success_json(
      project: @project.to_hash
    )
  end

  def projects
    projects = Project.all.map do |project|
      # Don't use proj.to_hash because we don't necessarily want to send back
      #   all the information.
      {
        project_name: project.project_name,
        project_name_full: project.project_name_full,
        resource: project.resource
      }
    end

    success_json({projects: projects})
  end

  def update_permission
    require_param(:email)
    @project = Project[project_name: @params[:project_name]]

    permission = @project.permissions.find do |p|
      p.user.email == @params[:email]
    end

    # fix strings from HTML POST
    @params[:privileged] = @params[:privileged] == 'true' if ['true', 'false'].include?(@params[:privileged])

    raise Etna::BadRequest, "No such user on project #{@params[:project_name]}!" unless permission

    raise Etna::Forbidden, 'Cannot update admin role!' if permission.role == 'administrator' && !@user.is_superuser?

    raise Etna::Forbidden, 'Cannot grant admin role!' if @params[:role] == 'administrator' && !@user.is_superuser?

    raise Etna::BadRequest, "Unknown role #{@params[:role]}" unless !@params[:role] || ['administrator', 'viewer', 'editor', 'disabled'].include?(@params[:role])

    if @params[:role] == 'disabled'
      permission.delete
    else
      permission.role = @params[:role] if @params[:role]
      permission.privileged = @params[:privileged] if [true, false].include?(@params[:privileged])
      permission.affiliation = @params[:affiliation] if @params[:affiliation]
      permission.save
    end

    @response.redirect("/#{@params[:project_name]}")
    @response.finish
  end

  def add_user
    require_params(:email, :name, :role)
    @project = Project[project_name: @params[:project_name]]

    @email = @params[:email].downcase.strip

    raise Etna::Forbidden, 'Cannot set admin role!' if @params[:role] == 'administrator'

    raise Etna::BadRequest, "Unknown role #{@params[:role]}" unless ['viewer', 'editor'].include?(@params[:role])

    unless @project.permissions.any? { |p| p.user.email == @email }
      user = User[email: @email]
      name = @params[:name]&.strip

      unless user
        raise Etna::BadRequest, 'Badly formed email address' unless @email =~ URI::MailTo::EMAIL_REGEXP

        raise Etna::BadRequest, 'Missing name' if name.empty?
        user = User.create(email: @email, name: name)
      end

      permission = Permission.create(project: @project, user: user, role: @params[:role])
      permission.role = @params[:role] if ['viewer', 'editor'].include?(@params[:role])
      permission.privileged = false
      permission.affiliation = @params[:affiliation]
      permission.save
    end

    @response.redirect("/#{@params[:project_name]}")
    @response.finish
  end

  def add_project
    require_params(:project_name, :project_name_full)

    raise Etna::BadRequest, "project_name should be like #{Project::PROJECT_NAME_MATCH.source}" unless Project.valid_name?(@params[:project_name])

    raise Etna::BadRequest, 'project_name_full cannot be empty' if @params[:project_name_full].nil? || @params[:project_name_full].empty?

    project = Project[project_name: @params[:project_name]]

    if project.nil?
      project = Project.create(
          project_name: @params[:project_name],
          project_name_full: @params[:project_name_full]
      )
    end

    @response.redirect('/')
    @response.finish
  end

  def update_project
    @project = Project[project_name: @params[:project_name]]

    raise Etna::BadRequest, "invalid project" if @project.nil?

    @project.update(resource: !!@params[:resource]) unless @params[:resource].nil?
    @project.refresh

    success_json(@project.to_hash)
  end

  def flag_user
    require_params(:flags, :email)

    user = User.find(email: @params[:email])

    if @params[:flags] &&
        !(@params[:flags].is_a?(Array) &&
            @params[:flags].all? { |f| f.is_a?(String) && f =~ /^\w+$/ })
      raise Etna::BadRequest, "Flags should be an array of words per /^\w+$/"
    end

    raise Etna::BadRequest, "No such user #{@params[:email]}" unless user

    user.update(flags: @params[:flags])

    success_json(user.to_hash)
  end
end
