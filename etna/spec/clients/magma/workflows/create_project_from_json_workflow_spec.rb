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
    stub_janus_setup
    # Magma retrieve needs to update the @target_models as we create them
    stub_request(:post, /#{MAGMA_HOST}\/retrieve/).
      to_return({body: {models: {}}.to_json}).then.
      to_return({body: {models: {
          assay_name: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp,
          project: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp,
          project: model_stamp,
          timepoint: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp,
          project: model_stamp,
          timepoint: model_stamp,
          patient: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp,
          project: model_stamp,
          timepoint: model_stamp,
          patient: model_stamp,
          document: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp,
          project: model_stamp,
          timepoint: model_stamp,
          patient: model_stamp,
          document: model_stamp,
          status: model_stamp
        }}.to_json}).times(3).then.
      to_return({body: {models: {
          assay_name: model_stamp,
          assay_pool: model_stamp,
          project: model_stamp,
          timepoint: model_stamp,
          patient: model_stamp,
          document: model_stamp,
          status: model_stamp,
          symptom: model_stamp
        }}.to_json})

    stub_magma_update_model
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
      with { |req| req.body.include?('add_attribute') }.times(60)
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
