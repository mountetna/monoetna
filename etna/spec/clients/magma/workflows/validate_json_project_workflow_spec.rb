require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::ValidateJsonProjectWorkflow do
  before(:each) do
  end

  it 'raises exception for invalid project JSON files' do
    workflow = Etna::Clients::Magma::ValidateJsonProjectWorkflow.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_missing_project_keys.json'
    )
    expect {
      workflow.validate
    }.to raise_error(Exception)
  end

  it 'does not raise exception for valid JSON files' do
    workflow = Etna::Clients::Magma::ValidateJsonProjectWorkflow.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )

    expect {
      workflow.validate
    }.not_to raise_error(Exception)
  end
end
