require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::AttributeActionsFromJsonWorkflow do
  let(:magma_client) {Etna::Clients::Magma.new(
    token: TEST_TOKEN,
    host: MAGMA_HOST)}

  before(:each) do
    @all_updates = []
    stub_magma_models(
      JSON.parse(File.read('./spec/fixtures/attribute_actions/test_project_magma_models.json')))
    stub_magma_update_model
  end

  def updates
    @all_updates.inject({}) do |acc, n|
      n.keys.each do |k|
        (acc[k] ||= {}).update(n[k])
      end
      acc
    end
  end

  it 'executes valid actions against Magma' do
    workflow = Etna::Clients::Magma::AttributeActionsFromJsonWorkflow.new(
      magma_client: magma_client,
      project_name: PROJECT,
      filepath: './spec/fixtures/attribute_actions/multiple_actions_valid.json'
    )
    workflow.run!

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: "Etna #{TEST_TOKEN}"}).
      with { |req| req.body.include?('add_attribute') }.times(1)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: "Etna #{TEST_TOKEN}"}).
      with { |req| req.body.include?('rename_attribute') }.times(1)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: "Etna #{TEST_TOKEN}"}).
      with { |req| req.body.include?('update_attribute') }.times(1)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: "Etna #{TEST_TOKEN}"}).
      with { |req| req.body.include?('add_link') }.times(1)
  end
end
