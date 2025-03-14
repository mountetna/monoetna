class AdminController < Janus::Controller
  def main
    @janus_user =  janus_user
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

    raise Etna::BadRequest, "Unknown role #{@params[:role]}" unless !@params[:role] || Token::ROLE_KEYS.values.concat(['disabled']).include?(@params[:role])

    if @params[:role] == 'disabled'
      permission.delete
    else
      permission.role = @params[:role] if @params[:role]
      permission.privileged = @params[:privileged] if [true, false].include?(@params[:privileged])
      permission.affiliation = @params[:affiliation] if @params[:affiliation]
      permission.save
    end

    event_log(
      event: 'update_permission',
      message: "changed #{@params[:email]} to #{@params.slice(:role,:affiliation,:privileged)}"
    )

    @response.redirect("/#{@params[:project_name]}")
    @response.finish
  end

  def add_user
    require_params(:email, :name, :role)

    settable_roles = ['viewer', 'editor', 'guest']

    @project = Project[project_name: @params[:project_name]]

    @email = @params[:email].downcase.strip

    raise Etna::Forbidden, 'Cannot set admin role!' if @params[:role] == 'administrator'

    raise Etna::BadRequest, "Unknown role #{@params[:role]}" unless settable_roles.include?(@params[:role])

    unless @project.permissions.any? { |p| p.user.email == @email }
      user = User[email: @email]
      name = @params[:name]&.strip

      unless user
        raise Etna::BadRequest, 'Badly formed email address' unless @email =~ URI::MailTo::EMAIL_REGEXP

        raise Etna::BadRequest, 'Missing name' if name.empty?
        user = User.create(email: @email, name: name)
      end

      permission = Permission.create(project: @project, user: user, role: @params[:role])
      permission.role = @params[:role] if settable_roles.include?(@params[:role])
      permission.privileged = false
      permission.affiliation = @params[:affiliation]
      permission.save

      event_log(
        event: 'add_user',
        message: "added user #{@email} with role #{@params[:role]}"
      )
    end

    @response.redirect("/#{@params[:project_name]}")
    @response.finish
  end

  def add_project
    require_params(:project_name, :project_name_full)

    raise Etna::BadRequest, "project_name should be like #{Project::PROJECT_NAME_MATCH.source}" unless Project.valid_name?(@params[:project_name])

    raise Etna::BadRequest, 'project_name_full cannot be empty' if @params[:project_name_full].nil? || @params[:project_name_full].empty?

    raise Etna::BadRequest, 'template_project_name does not exist' if template_project_name_provided? && !template_project_exists?

    project = Project[project_name: @params[:project_name]]

    if project.nil?
      project = Project.create(
          project_name: @params[:project_name],
          project_name_full: @params[:project_name_full]
      )

      Permission.create(project: project, user: janus_user, role: 'administrator')

      janus_user.refresh
      refreshed_token = janus_user.create_token!

      magma_client = magma_client(refreshed_token)

      begin
        magma_client.update_model(
          Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project.project_name,
            actions: [
              Etna::Clients::Magma::AddProjectAction.new
            ]
          )
        )

        magma_client.update_json(Etna::Clients::Magma::UpdateRequest.new(
          project_name: project.project_name,
          revisions: {
              'project' => { project.project_name => { name: project.project_name } },
          }))

        copy_project_template(magma_client) if template_project_name_provided?
      end if magma_client

      event_log(
        event: 'create_project',
        message: "created #{@params[:project_name_full]}"
      )
      Janus.instance.set_token_cookie(@response, refreshed_token)
    end

    @response.redirect('/')
    @response.finish
  end

  def update_project
    @project = Project[project_name: @params[:project_name]]

    raise Etna::BadRequest, "invalid project" if @project.nil?
    raise Etna::BadRequest, "invalid contact email" if @params[:contact_email] && !valid_contact?

    @project.update(**update_payload)

    event_log(
      event: 'update_project',
      message: "updated project settings",
      payload: update_payload
    )

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

  def update_cc_agreement
    # User agrees to the code of conduct for the specified project
    require_params(:project_name, :cc_text, :agreed)
    project_name = @params[:project_name]
    project = Project[project_name: @params[:project_name]]
    if project.nil?
      return failure(404, "Project #{project_name} does not exist.")
    end

    if !project.requires_agreement
      return failure(403, "Project #{project_name} does not require a code of conduct agreement.")
    end

    agreement = CcAgreement.create(
      user_email: @user.email,
      project_name: project_name,
      cc_text: @params[:cc_text],
      agreed: !!@params[:agreed]
    )

    janus_user.set_guest_permissions!(project_name)

    event_log(
      event: 'cc_agreement',
      message: "signed guest code of conduct agreement"
    )

    if project.contact_email && project.contact_email.length > 0 && @params[:agreed]
      send_email(
        "Project Lead of #{project_name}",
        project.contact_email,
        "#{project_name} Code of Conduct was signed by #{@user.email}",
        "User #{@user.email} has agreed to the Code of Conduct and will be enjoying guest level access to the project."
      )
    end

    token = janus_user.create_token!
    Janus.instance.set_token_cookie(@response, token)
    success_json(agreement.to_hash)
  end

  private

  def update_payload
    {}.tap do |result|
      result[:resource] = !!@params[:resource] unless @params[:resource].nil?
      result[:requires_agreement] = !!@params[:requires_agreement] unless @params[:requires_agreement].nil?
      result[:cc_text] = @params[:cc_text] unless @params[:cc_text].nil?
      result[:contact_email] = @params[:contact_email].strip if valid_contact?
      result[:project_name_full] = @params[:project_name_full] if @params[:project_name_full]&.present?
    end
  end

  def valid_contact?
    @params[:contact_email]&.empty? || @params[:contact_email]&.strip&.rpartition('@')&.last == Janus.instance.config(:token_domain)
  end

  def magma_client(token)
    return nil unless Janus.instance.config(:magma) && Janus.instance.config(:magma)[:host]

    Etna::Clients::Magma.new(
      host: Janus.instance.config(:magma)[:host],
      token: token
    )
  end

  def janus_user
    User[email: @user.email]
  end

  def template_project_name_provided?
    !!@params[:template_project_name]
  end

  def template_project_exists?
    return false unless template_project_name_provided?

    !Project[project_name: @params[:template_project_name]].nil?
  end

  def copy_project_template(target_client)
    # Do this in a separate thread ... can take more than 1 minute, so
    #   the connection times out.
    Thread.new do
      begin
        source_models = target_client.retrieve(
          Etna::Clients::Magma::RetrievalRequest.new(
            project_name: @params[:template_project_name]
          )
        ).models

        workflow = Etna::Clients::Magma::ModelSynchronizationWorkflow.new(
          target_client: target_client,
          source_models: source_models,
          target_project: @params[:project_name]
        )
        workflow.ensure_model_tree('project')
      rescue => e
        Janus.instance.logger.error("Error in thread for copying project template! #{e.message}")
        Janus.instance.logger.log_error(e)
      end
    end
  end
end
