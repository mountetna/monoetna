describe AdminController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#projects' do
    before(:each) do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror', resource: true)
    end

    it 'prevents access to the list of all projects for non-superusers' do
      # ordinary user cannot
      auth_header(:janus)
      get('/allprojects')

      expect(last_response.status).to eq(403)
    end

    it 'returns a list of all projects' do
      auth_header(:superuser)
      get('/allprojects')
      expect(last_response.status).to eq(200)

      expect(json_body[:projects]).to eq(
        [
          { project_name: "gateway", project_name_full: "Gateway", resource: false},
          { project_name: "tunnel", project_name_full: "Tunnel", resource: false},
          { project_name: "mirror", project_name_full: "Mirror", resource: true}
        ]
      )
    end
  end

  context '#project' do
    it 'returns a project view to the admin' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator')

      auth_header(:janus)
      get('/door')

      expect(last_response.status).to eq(200)
    end

    it 'forbids the project view to viewers' do
      user = create(:user, name: 'Lar Familiaris', email: 'lar@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'viewer')

      auth_header(:lar)
      get('/door')

      expect(last_response.status).to eq(403)
    end

    it 'forbids project outsiders' do
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      auth_header(:viewer)
      get('/door')

      expect(last_response.status).to eq(403)
    end

    it 'viewer cannot set the resource flag' do
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      auth_header(:viewer)
      json_post('door', resource: true)

      expect(last_response.status).to eq(403)
    end

    it 'can set the resource flag' do
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      expect(door.resource).to eq(false)

      auth_header(:portunus)

      json_post('door', resource: true)

      expect(last_response.status).to eq(200)
      door.refresh
      expect(door.resource).to eq(true)

      json_post('door', resource: false)

      expect(last_response.status).to eq(200)
      door.refresh
      expect(door.resource).to eq(false)
    end

    it 'resource flag does not change if not provided in params' do
      door = create(:project, project_name: 'door', project_name_full: 'Door', resource: true)

      expect(door.resource).to eq(true)

      auth_header(:portunus)
      json_post('door', arbitrary: false)

      expect(last_response.status).to eq(200)
      door.refresh
      expect(door.resource).to eq(true)
    end

    it 'returns a list of permissions for the project' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      get('/project/door')

      expect(last_response.status).to eq(200)
      expect(json_body[:project]).to match(
        permissions: [
          {
            affiliation: nil,
            privileged: true,
            project_name: "door",
            role: "administrator",
            user_email: "janus@two-faces.org",
            user_name: "Janus Bifrons"
          },
          {
            affiliation: nil,
            privileged: nil,
            project_name: "door",
            role: "editor",
            user_email: "portunus@two-faces.org",
            user_name: "Portunus"
          }
        ],
        project_name: "door",
        project_name_full: "Door",
        resource: false
      )
    end
  end

  context '#update_permission' do
    it 'allows an admin to update a role' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'viewer')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      perm2.refresh
      expect(perm2.role).to eq('viewer')
    end

    it 'allows an admin to update an affiliation' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', affiliation: 'ILWU')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      perm2.refresh
      expect(perm2.affiliation).to eq('ILWU')
    end

    it 'allows an admin to grant privilege' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', privileged: true)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      perm2.refresh
      expect(perm2).to be_privileged
    end

    it 'allows an admin to remove privilege' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor', privileged: true)

      auth_header(:janus)
      json_post('update_permission/door', email: 'portunus@two-faces.org', privileged: false)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      perm2.refresh
      expect(perm2).not_to be_privileged
    end

    it 'forbids a non-admin from updating roles' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'editor')

      auth_header(:superuser)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'administrator')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      perm2.refresh
      expect(perm2.role).to eq('administrator')
    end

    it 'allows the superuser to remove admin powers' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'administrator')

      auth_header(:superuser)
      json_post('update_permission/door', email: 'portunus@two-faces.org', role: 'editor')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      perm2.refresh
      expect(perm2.role).to eq('editor')
    end

    it 'deletes a permission' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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

  context '#add_user' do
    it 'allows an admin to add a new user to a project' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor', affiliation: "ILWU")

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

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

    it 'strips leading and trailing whitespace from an e-mail' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: ' portunus@two-faces.org ', name: 'Portunus', role: 'editor', affiliation: "ILWU")

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor')

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

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

    it 'allows an admin to add a user to a project twice' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)
      perm2 = create(:permission, project: door, user: user2, role: 'viewer', privileged: false)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor', privileged: true)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

      expect(Permission.count).to eq(2)
      expect(User.count).to eq(2)

      user2.refresh
      perm2.refresh

      # nothing has changed
      expect(perm2).not_to be_privileged
      expect(perm2.user).to eq(user2)
      expect(perm2.project).to eq(door)
      expect(perm2.role).to eq('viewer')

      expect(user2.name).to eq('Portunus')
      expect(user2.email).to eq('portunus@two-faces.org')
    end

    it 'does not allow admin to give privilege to a new user' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'portunus@two-faces.org', name: 'Portunus', role: 'editor', privileged: true)

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq('/door')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

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

    it 'squashes case in email addresses' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      door = create(:project, project_name: 'door', project_name_full: 'Door')
      perm = create(:permission, project: door, user: user, role: 'administrator', privileged: true)

      auth_header(:janus)
      json_post('add_user/door', email: 'Portunus@two-faces.org', name: 'Portunus', role: 'editor', privileged: true)

      expect(last_response.status).to eq(302)

      expect(Permission.count).to eq(2)
      expect(User.count).to eq(2)
      expect(User.first.email).to eq(user.email)
      expect(User.last.email).to eq('portunus@two-faces.org')
    end

    it 'forbids a non-admin from adding a user' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      #user2 = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

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

  context '#add_project' do
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
      [ 'Doors', "Doo\nrs", ' doors', 'doors	' , '1x_door', 'pg_door', 'door_2_project'].each do |name|
        json_post('add_project', project_name: name, project_name_full: name)

        expect(last_response.status).to eq(422)
        expect(json_body[:error]).to match(/project_name should be like/)

        expect(Project.count).to eq(0)
      end
    end

    it 'requires some project_name_full' do
      auth_header(:superuser)
      json_post('add_project', project_name: 'door', project_name_full: '')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('project_name_full cannot be empty')

      expect(Project.count).to eq(0)
    end
  end

  context '#flag_user' do
    it 'sets flags on the user' do
      user = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      auth_header(:superuser)
      json_post('flag_user', email: 'portunus@two-faces.org', flags: [ 'doors' ])

      # the flags are set
      expect(last_response.status).to eq(200)
      expect(json_body[:flags]).to eq([ 'doors' ])

      # we get the flags back
      user.refresh
      expect(user.flags).to eq([ 'doors' ])
    end

    it 'clears flags on the user' do
      user = create(:user, name: 'Portunus', email: 'portunus@two-faces.org', flags: [ 'doors' ])

      auth_header(:superuser)
      json_post('flag_user', email: 'portunus@two-faces.org', flags: nil)

      # the flags are set
      expect(last_response.status).to eq(200)
      expect(json_body[:flags]).to eq(nil)

      # we get the flags back
      user.refresh
      expect(user.flags).to eq(nil)
    end

    it 'does not set invalid flags' do
      user = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      auth_header(:superuser)
      json_post('flag_user', email: 'portunus@two-faces.org', flags: [ 'lll', 2 ])

      expect(last_response.status).to eq(422)

      # the flags are unchanged
      user.refresh
      expect(user.flags).to eq(nil)
    end

    it 'prevents non-superusers from setting flags' do
      user = create(:user, name: 'Portunus', email: 'portunus@two-faces.org')

      auth_header(:portunus)
      json_post('flag_user', email: 'portunus@two-faces.org', flags: [ 'doors' ])

      expect(last_response.status).to eq(403)

      # the flags are not changed
      user.refresh
      expect(user.flags).to eq(nil)
    end
  end
end
