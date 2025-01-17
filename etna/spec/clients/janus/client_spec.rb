require 'webmock/rspec'
require 'json'
require_relative '../../../lib/etna/clients/janus'

describe 'Janus Client class' do
  let(:test_class) { Etna::Clients::Janus.new(token: TEST_TOKEN, host: JANUS_HOST) }

  before(:each) do
    stub_janus_setup
  end

  it 'can fetch project details' do
    test_class.get_project(Etna::Clients::Janus::GetProjectRequest.new(
      project_name: 'test'
    ))
    expect(WebMock).to have_requested(:get, %r!#{JANUS_HOST}/api/admin/test!)
  end

  it 'can fetch user projects' do
    test_class.get_projects()
    expect(WebMock).to have_requested(:get, %r!#{JANUS_HOST}/api/user/projects!)
  end

  it 'can add a new project' do
    test_class.add_project(Etna::Clients::Janus::AddProjectRequest.new(
      project_name: 'test',
      project_name_full: 'TestProject1'
    ))
    expect(WebMock).to have_requested(:post, %r!#{JANUS_HOST}/api/admin/add_project!)
  end

  it 'can add a new user to a project' do
    test_class.add_user(Etna::Clients::Janus::AddUserRequest.new(
      project_name: 'test',
      email: 'tester@janus.test',
      role: 'viewer',
      affiliation: 'None'
    ))

    expect(WebMock).to have_requested(:post, %r!#{JANUS_HOST}/api/admin/test/add_user!)
  end

  it 'can update permission for a project' do
    test_class.update_permission(Etna::Clients::Janus::UpdatePermissionRequest.new(
      project_name: 'test',
      email: 'tester@janus.test',
      role: 'viewer',
      affiliation: 'None'
    ))

    expect(WebMock).to have_requested(:post, %r!#{JANUS_HOST}/api/admin/test/update_permission!)
  end

  it 'can refresh a user\'s token' do
    response = test_class.refresh_token(Etna::Clients::Janus::RefreshTokenRequest.new)
    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/api\/tokens\/generate/)
    expect(response.token).to eq('a token for you!')
  end

  it 'can fetch stats' do
    stub_janus_get_project_stats
    test_class.get_project_stats(Etna::Clients::Janus::GetStatsRequest.new)
    expect(WebMock).to have_requested(:get, %r!#{JANUS_HOST}/api/stats/projects!)
  end

  it 'can fetch stats with projects specified' do
    project_names = ['test']

    stub_janus_get_project_stats(project_names)

    query = project_names.map do |name|
      "projects%5B%5D=#{name}"
    end.join('&')

    test_class.get_project_stats(Etna::Clients::Janus::GetStatsRequest.new(project_names: project_names))
    expect(WebMock).to have_requested(:get, %r!#{JANUS_HOST}/api/stats/projects\?#{query}!)
  end
end
