require 'webmock/rspec'
require 'json'
require_relative '../../../../lib/etna/clients/janus'
require 'pry'

def model_stamp
  {
    template: {
      attributes: {}
    }
  }
end

describe Etna::Clients::Magma::CreateProjectFromJsonWorkflow do
  let(:magma_client) {Etna::Clients::Magma.new(
    token: '123',
    host: MAGMA_HOST)}
  let(:janus_client) {Etna::Clients::Janus.new(
    token: '123',
    host: JANUS_HOST)}

  before(:each) do
    @all_updates = []
    stub_janus_setup
    stub_magma_models({
      models: {
      assay_name: model_stamp,
      assay_pool: model_stamp,
      project: model_stamp,
      timepoint: model_stamp,
      patient: model_stamp,
      document: model_stamp,
      status: model_stamp,
      symptom: model_stamp
    }})
    stub_magma_update_model
    stub_magma_update
  end

  def updates
    @all_updates.inject({}) do |acc, n|
      n.keys.each do |k|
        (acc[k] ||= {}).update(n[k])
      end
      acc
    end
  end

  it 'creates the project in janus and magma' do
    workflow = Etna::Clients::Magma::CreateProjectFromJsonWorkflow.new(
      magma_client: magma_client,
      janus_client: janus_client,
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )
    workflow.create!

    expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/add_project/).
      with(headers: {Authorization: 'Etna 123'})
      expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/add_user\/#{PROJECT}/)
      expect(WebMock).to have_requested(:post, /#{JANUS_HOST}\/update_permission\/#{PROJECT}/)
    expect(WebMock).to have_requested(:get, /#{JANUS_HOST}\/refresh_token/)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna a token for you!'}).
      with { |req| req.body.include?('add_model') }.times(8)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna a token for you!'}).
      with { |req| req.body.include?('add_attribute') }.times(62)
    # The call above gets executed a lot of times, because our Mock of the
    #   /retrieve endpoint is static. So it doesn't reflect newly created
    #   attributes. For the models higher in the tree (closer to the root),
    #   this means that ensure_model_tree attempts to add their
    #   attributes multiple times. This does not happen in a non-test
    #   environment.

    # Make sure the assay_name identifier validation is submitted.
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna a token for you!'}).
      with { |req| req.body.include?('update_attribute') }

    expect(updates).to eq({
      "project" => {
        "test" => {"name" => "test"}
      }
    })
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
