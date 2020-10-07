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

  it 'reports errors for invalid attribute information' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_attribute_errors.json'
    )

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(16)
    expect(project.errors).to eq([
      "Missing required key for model assay_name, attribute tube_name, validation: \"value\".",
      "Missing required key for model assay_name, attribute assay_pool: \"desc\".",
      "Missing required key for model assay_name, attribute vendor: \"desc\".",
      "Invalid type for model assay_name, attribute vendor, validation: \"Lo que sea\".",
      "Invalid attribute_type for model document, attribute document_desc: \"\".",
      "Invalid attribute_type for model document, attribute document_desc: \"\".",
      "Invalid attribute_type for model document, attribute version: \"copacetic\".",
      "Invalid attribute_type for model document, attribute version_date: \"date\".",
      "Missing required key for model assay_pool, attribute biospecimen: \"attribute_type\".",
      "Missing required key for model assay_pool, attribute biospecimen: \"desc\".",
      "Invalid value for model assay_pool, attribute biospecimen, validation: \"\".",
      "Missing required key for model assay_pool, attribute cells_loaded: \"display_name\".",
      "Missing required key for model assay_pool, attribute project: \"desc\".",
      "Missing required key for model timepoint, attribute day: \"desc\".",
      "Linked model, \"pool_deep_end\", on attribute assay_pool of model assay_name does not exist!",
      "Model \"assay_pool\" already belongs to parent model \"project\". Remove attribute \"project\"."
    ])
  end

  it 'can get models as Magma::Models' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/create_project_fixture_valid.json'
    )

    source_models = project.get_magma_models

    # Make sure the identifier, parent, and linking attributes were translated correctly!
    project = source_models.model('project')
    expect(project.template.attributes.attribute('document').attribute_type).to eq(Etna::Clients::Magma::AttributeType::COLLECTION)
    expect(project.template.attributes.attribute('patient').attribute_type).to eq(Etna::Clients::Magma::AttributeType::COLLECTION)
    expect(project.template.attributes.attribute('assay_pool').attribute_type).to eq(Etna::Clients::Magma::AttributeType::COLLECTION)

    assay_pool = source_models.model('assay_pool')
    expect(assay_pool.template.attributes.attribute('assay_name').attribute_type).to eq(Etna::Clients::Magma::AttributeType::COLLECTION)
    expect(assay_pool.template.attributes.attribute('assay_name').link_model_name).to eq('assay_name')
    expect(assay_pool.template.attributes.attribute('project').attribute_type).to eq(Etna::Clients::Magma::AttributeType::PARENT)

    assay_name = source_models.model('assay_name')
    expect(assay_name.template.attributes.attribute('assay_pool').attribute_type).to eq(Etna::Clients::Magma::AttributeType::LINK)
    expect(assay_name.template.attributes.attribute('timepoint').attribute_type).to eq(Etna::Clients::Magma::AttributeType::PARENT)
    expect(assay_name.template.attributes.attribute('tube_name').attribute_type).to eq(Etna::Clients::Magma::AttributeType::IDENTIFIER)

    status = source_models.model('status')
    expect(status.template.attributes.attribute('patient').attribute_type).to eq(Etna::Clients::Magma::AttributeType::PARENT)
    expect(status.template.attributes.attribute('patient').link_model_name).to eq('patient')
  end
end
