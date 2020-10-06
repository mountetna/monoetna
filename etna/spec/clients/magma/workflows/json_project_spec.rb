require 'webmock/rspec'
require 'json'
require 'pry'

describe Etna::Clients::Magma::JsonProject do
  before(:each) do
  end

  it 'raises exception for missing project keys' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_missing_project_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(2)
    expect(project.errors.first).to eq('Missing required key for root project, "project_name".')
    expect(project.errors.last).to eq('Missing required key for root project, "project_name_full".')
  end

  it 'loads a project JSON definition' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )

    expect(project.valid?).to eq(true)
    expect(project.name).to eq('test1')
    expect(project.name_full).to eq('Testing your immune system')
  end
end
