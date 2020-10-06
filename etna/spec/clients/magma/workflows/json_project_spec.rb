require 'webmock/rspec'
require 'json'
require 'pry'

describe Etna::Clients::Magma::JsonProject do
  before(:each) do
  end

  it 'raises exception for missing project name' do
    expect {
      Etna::Clients::Magma::JsonProject.new(
        filepath: './spec/fixtures/create_project/create_project_fixture_missing_project_name.json'
      )
    }.to raise_error(
      RuntimeError,
      'Missing required key, "project_name".')
  end

  it 'raises exception for missing project_name_full' do
    expect {
      Etna::Clients::Magma::JsonProject.new(
        filepath: './spec/fixtures/create_project/create_project_fixture_missing_project_name_full.json'
      )
    }.to raise_error(
      RuntimeError,
      'Missing required key, "project_name_full".')
  end

  it 'loads a project JSON definition' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )

    expect(project.name).to eq('test1')
    expect(project.name_full).to eq('Testing your immune system')
  end
end
