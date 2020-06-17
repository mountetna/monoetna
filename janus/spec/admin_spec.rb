describe AdminController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context 'main' do
    it 'returns a list of user projects' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')

      perm = create(:permission, project: mirror, user: user, role: 'editor')
      perm = create(:permission, project: gateway, user: user, role: 'editor')

      auth_header(:janus)
      get('/')

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/Your Projects/)
      expect(last_response.body).to match(/Gateway/)
      expect(last_response.body).to match(/Mirror/)
      expect(last_response.body).not_to match(/Tunnel/)
    end

    it 'returns the user public key fingerprint' do
      pkey = OpenSSL::PKey::RSA.new(1024)
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org', public_key: pkey.public_key)

      auth_header(:janus)

      get('/')

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/Your Keys/)
      expect(last_response.body).to match(/#{user.key_fingerprint}/i)
    end
  end

  context 'project' do
    it 'returns a project view to the admin' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator')

      auth_header(:janus)
      get('/project/door')

      expect(last_response.status).to eq(200)
    end

    it 'shows a static project view to editors' do
      user = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'editor')

      auth_header(:portunus)
      get('/project/door')

      expect(last_response.status).to eq(200)
      expect(html_body.css('input')).to be_empty
    end

    it 'forbids the project view to viewers' do
      user = create(:user, first_name: 'Lar', last_name: 'Familiaris', email: 'lar@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'viewer')

      auth_header(:lar)
      get('/project/door')

      expect(last_response.status).to eq(403)
    end

    it 'forbids project outsiders' do
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      auth_header(:viewer)
      get('/project/door')

      expect(last_response.status).to eq(403)
    end

    it 'returns a list of permissions for the project' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      get('/project/door')

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/Door/)
    end
  end

  context 'update_permission' do
    it 'allows an admin to update a role' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'viewer')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      perm2.refresh
      expect(perm2.role).to eq('viewer')
    end

    it 'allows an admin to update an affiliation' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', affiliation: 'ILWU')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      perm2.refresh
      expect(perm2.affiliation).to eq('ILWU')
    end

    it 'allows an admin to grant privilege' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', privileged: true)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      perm2.refresh
      expect(perm2).to be_privileged
    end

    it 'allows an admin to remove privilege' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor', privileged: true)

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', privileged: false)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      perm2.refresh
      expect(perm2).not_to be_privileged
    end

    it 'forbids a non-admin from updating roles' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:portunus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'viewer')

      expect(last_response.status).to eq(403)

      perm2.refresh
      expect(perm2.role).to eq('editor')
    end

    it 'forbids a non-admin from updating affiliation' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:portunus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', affiliation: 'ILWU')

      expect(last_response.status).to eq(403)

      perm2.refresh
      expect(perm2.affiliation).to be_nil
    end

    it 'forbids a non-admin from updating privileges' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:portunus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', privileged: true)

      expect(last_response.status).to eq(403)

      perm2.refresh
      expect(perm2).not_to be_privileged
    end

    it 'forbids an admin from updating role for an admin' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'administrator')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'editor')

      expect(last_response.status).to eq(403)

      perm2.refresh
      expect(perm2.role).to eq('administrator')
    end

    it 'allows the superuser to grant admin powers' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:superuser)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'administrator')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      perm2.refresh
      expect(perm2.role).to eq('administrator')
    end

    it 'allows the superuser to remove admin powers' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'administrator')

      auth_header(:superuser)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'editor')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      perm2.refresh
      expect(perm2.role).to eq('editor')
    end

    it 'deletes a permission' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'disabled')

      expect(last_response.status).to eq(302)

      door.refresh
      expect(door.permissions).to eq([perm])
    end
  end

  context 'add_user' do
    it 'allows an admin to add a new user to a project' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor', affiliation: "ILWU")

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      expect(Permission.count).to eq(2)
      expect(User.count).to eq(2)

      perm2 = Permission.last
      user2 = User.last
      expect(perm2).not_to be_privileged
      expect(perm2.user).to eq(user2)
      expect(perm2.project).to eq(door)
      expect(perm2.role).to eq('editor')
      expect(perm2.affiliation).to eq('ILWU')

      expect(user2.name).to eq('Portunus')
      expect(user2.email).to eq('portunus@two-faces.org')
    end

    it 'allows an admin to add an existing user to a project' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      expect(Permission.count).to eq(2)
      expect(User.count).to eq(2)

      perm2 = Permission.last
      user2 = User.last
      expect(perm2).not_to be_privileged
      expect(perm2.user).to eq(user2)
      expect(perm2.project).to eq(door)
      expect(perm2.role).to eq('editor')

      expect(user2.name).to eq('Portunus')
      expect(user2.email).to eq('portunus@two-faces.org')
    end

    it 'does not allow admin to give privilege to a new user' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor', privileged: true)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/project/door')

      expect(Permission.count).to eq(2)
      expect(User.count).to eq(2)

      perm2 = Permission.last
      user2 = User.last
      expect(perm2).not_to be_privileged
      expect(perm2.user).to eq(user2)
      expect(perm2.project).to eq(door)
      expect(perm2.role).to eq('editor')

      expect(user2.name).to eq('Portunus')
      expect(user2.email).to eq('portunus@two-faces.org')
    end

    it 'does not allow admin to give admin permission' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'administrator')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot set admin role!')

      expect(Permission.count).to eq(1)
      expect(User.count).to eq(2)
      expect(Permission.first.user).to eq(user)
    end

    it 'rejects incorrect email addresses for new users' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus.two-faces.org', name: 'Portunus', role: 'editor', privileged: true)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Badly formed email address')

      expect(Permission.count).to eq(1)
      expect(User.count).to eq(1)
      expect(User.first).to eq(user)
      expect(Permission.first.user).to eq(user)
    end

    it 'forbids a non-admin from adding a user' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      #user2 = create(:user, first_name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      #perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:portunus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor')

      expect(last_response.status).to eq(403)

      expect(Permission.count).to eq(1)
      expect(User.count).to eq(1)
      expect(User.first).to eq(user)
      expect(Permission.first.user).to eq(user)
    end
  end

  context 'add_project' do
    it 'allows a superuser to add a new project' do
      auth_header(:superuser)
      json_post('add_project', project_name: 'door', project_name_full: "Doors")

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/')

      expect(Project.count).to eq(1)
      project = Project.first

      expect(project.project_name).to eq('door')
    end

    it 'does not allow admins to add a new project' do
      auth_header(:janus)
      json_post('add_project', project_name: 'door', project_name_full: "Doors")

      expect(last_response.status).to eq(403)

      expect(Project.count).to eq(0)
    end

    it 'requires a well-formed project_name' do
      auth_header(:superuser)
      json_post('add_project', project_name: 'Door', project_name_full: "Doors")

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to match(/project_name should be like/)

      expect(Project.count).to eq(0)
    end

    it 'requires some project_name_full' do
      auth_header(:superuser)
      json_post('add_project', project_name: 'door', project_name_full: '')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('project_name_full cannot be empty')

      expect(Project.count).to eq(0)
    end
  end
end
