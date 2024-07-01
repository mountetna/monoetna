describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#project_stats" do
    it "returns stats for all projects if none specified" do
      user1 = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      user2 = create(:user, name: 'Vesta Bule', email: 'vesta@two-faces.org')

      admin = create(:project, project_name: 'administration', project_name_full: 'Administration')
      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')
      door = create(:project, project_name: 'door', project_name_full: 'Door', resource: true)

      create(:permission, project: admin, user: user1, role: 'editor')
      create(:permission, project: tunnel, user: user1, role: 'viewer')
      create(:permission, project: tunnel, user: user2, role: 'viewer')
      create(:permission, project: mirror, user: user1, role: 'editor')
      create(:permission, project: gateway, user: user1, role: 'editor')

      header('Authorization', "Etna #{user1.create_token!}")
      get('/api/stats/projects')

      expect(last_response.status).to eq(200)
      expect(json_body[:projects]).to eq([
        {project_name: "administration", project_name_full: "Administration", resource: false, user_count: 1},
        {project_name: "gateway", project_name_full: "Gateway", resource: false, user_count: 1},
        {project_name: "tunnel", project_name_full: "Tunnel", resource: false, user_count: 2},
        {project_name: "mirror", project_name_full: "Mirror", resource: false, user_count: 1},
        {project_name: "door", project_name_full: "Door", resource: true, user_count: 0},
      ])
      expect(json_body[:user_count]).to eq(2)
    end

    it "returns stats only for specified projects" do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      admin = create(:project, project_name: 'administration', project_name_full: 'Administration')
      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      create(:permission, project: admin, user: user, role: 'editor')
      create(:permission, project: tunnel, user: user, role: 'viewer')
      create(:permission, project: mirror, user: user, role: 'editor')
      create(:permission, project: gateway, user: user, role: 'editor')


      header('Authorization', "Etna #{user.create_token!}")
      get('/api/stats/projects', projects: [tunnel.project_name, door.project_name])

      expect(last_response.status).to eq(200)
      expect(json_body[:projects]).to eq([
        {project_name: "tunnel", project_name_full: "Tunnel", resource: false, user_count: 1},
        {project_name: "door", project_name_full: "Door", resource: false, user_count: 0}
      ])
      expect(json_body[:user_count]).to eq(1)
    end

    it "only allows supereditors" do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')

      create(:permission, project: gateway, user: user, role: 'editor')

      header('Authorization', "Etna #{user.create_token!}")
      get('/api/stats/projects')

      expect(last_response.status).to eq(403)
    end
  end

end
