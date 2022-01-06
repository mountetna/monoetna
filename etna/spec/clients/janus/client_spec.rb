require 'webmock/rspec'
require 'json'
require_relative '../../../lib/etna/clients/janus'

describe 'Janus Client class' do
  let(:test_class) { Etna::Clients::Janus.new(token: TEST_TOKEN, host: 'https://janus.test') }

  before(:each) do
    stub_janus_setup
  end

  it 'can fetch project details' do
    test_class.get_project(Etna::Clients::Janus::GetProjectRequest.new(
      project_name: 'test'
    ))
    expect(WebMock).to have_requested(:get, /#{JANUS_HOST}\/project\/test/)
  end

  it 'can add a new project' do
    test_class.add_project(Etna::Clients::Janus::AddProjectRequest.new(
      project_name: 'test',
      project_name_full: 'TestProject1'
    ))
    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/add_project/)
  end

  it 'can add a new user to a project' do
    test_class.add_user(Etna::Clients::Janus::AddUserRequest.new(
      project_name: 'test',
      email: 'tester@janus.test',
      role: 'viewer',
      affiliation: 'None'
    ))

    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/add_user\/test/)
  end

  it 'can update permission for a project' do
    test_class.update_permission(Etna::Clients::Janus::UpdatePermissionRequest.new(
      project_name: 'test',
      email: 'tester@janus.test',
      role: 'viewer',
      affiliation: 'None'
    ))

    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/update_permission\/test/)
  end

  it 'can refresh a user\'s token' do
    response = test_class.refresh_token(Etna::Clients::Janus::RefreshTokenRequest.new)
    expect(WebMock).to have_requested(:get, /#{JANUS_HOST}\/refresh_token/)
    expect(response.token).to eq('a token for you!')
  end

  it 'can fetch a viewer-only version of a user\'s token' do
    response = test_class.viewer_token(Etna::Clients::Janus::ViewerTokenRequest.new)
    expect(WebMock).to have_requested(:get, /#{JANUS_HOST}\/viewer_token/)
    expect(response.token).to eq('a view-only token for you!')
  end
end