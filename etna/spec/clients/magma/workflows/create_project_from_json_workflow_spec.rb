require 'webmock/rspec'
require 'json'
require_relative '../../../../lib/etna/clients/janus'
require 'pry'

describe Etna::Clients::Magma::CreateProjectFromJsonWorkflow do
  let(:magma_client) {Etna::Clients::Magma.new(
    token: '123',
    host: MAGMA_HOST)}
  let(:janus_client) {Etna::Clients::Janus.new(
    token: '123',
    host: JANUS_HOST)}

  before(:each) do
    stub_janus_setup
    stub_magma_update_model
  end

  it 'creates the project in janus and magma' do
    workflow = Etna::Clients::Magma::CreateProjectFromJsonWorkflow.new(
      magma_client: magma_client,
      janus_client: janus_client,
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )
    workflow.create!

    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/add_project/)
    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/update_permission\/#{PROJECT}/)
    expect(WebMock).to have_requested(:get, /#{JANUS_HOST}\/refresh_token/)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).times(17) # 9 models, 8 with attributes
  end

  it 'raises exception for an invalid project JSON' do
    expect {
      Etna::Clients::Magma::CreateProjectFromJsonWorkflow.new(
        magma_client: magma_client,
        janus_client: janus_client,
        filepath: './spec/fixtures/create_project/create_project_fixture_missing_project_keys.json'
      )
    }.to raise_error(Exception)
  end
end
