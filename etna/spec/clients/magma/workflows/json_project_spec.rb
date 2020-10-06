require 'webmock/rspec'
require 'json'
require 'pry'

describe Etna::Clients::Magma::JsonProject do
  before(:each) do
  end

  it 'reports errors for missing project keys' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_missing_project_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(2)
    expect(project.errors.first).to eq('Missing required key for root project: "project_name".')
    expect(project.errors.last).to eq('Missing required key for root project: "project_name_full".')
  end

  it 'reports errors for blank project keys' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_blank_project_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(2)
    expect(project.errors.first).to eq('Invalid project_name for root project: "".')
    expect(project.errors.last).to eq('Invalid project_name_full for root project: "".')
  end

  it 'loads a project JSON definition' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )

    expect(project.valid?).to eq(true)
    expect(project.name).to eq('test1')
    expect(project.name_full).to eq('Testing your immune system')
  end

  it 'reports errors for missing model information' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_missing_model_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(6)
    expect(project.errors).to eq([
      "Missing required key for model assay_name: \"parent_link_type\".",
      "Missing required key for model document: \"parent_model_name\".",
      "Missing required key for model document: \"parent_link_type\".",
      "Missing required key for model assay_pool: \"identifier\".",
      "Missing required key for model patient: \"parent_model_name\".",
      "Missing required key for model project: \"identifier\"."
    ])
  end

  it 'reports errors for blank model information' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_blank_model_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(6)
    expect(project.errors).to eq([
      "Invalid parent_model_name for model assay_name: \"\".",
      "Invalid identifier for model document: \"\".",
      "Invalid parent_link_type for model assay_pool: \"\".",
      "Invalid parent_link_type for model assay_pool: \"\".",
      "Invalid parent_link_type for model patient: \"thingamajig\".",
      "Invalid identifier for model project: \"\"."
    ])
  end
end
