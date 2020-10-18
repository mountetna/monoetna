require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::JsonProject do
  before(:each) do
  end

  it 'reports errors for missing project keys' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/missing_project_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(3)
    expect(project.errors).to eq([
      'Missing required key for root project: "project_name".',
      'Missing required key for root project: "project_name_full".',
      'Project name  must be snake_case and cannot start with a number or "pg_".'
    ])
  end

  it 'reports errors for blank project keys' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/blank_project_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(3)
    expect(project.errors).to eq([
      'Invalid empty project_name for root project: "".',
      'Invalid empty project_name_full for root project: "".',
      'Project name  must be snake_case and cannot start with a number or "pg_".'
    ])
  end

  it 'reports errors for invalid project name' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/invalid_project_name.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(1)
    expect(project.errors.first).to eq('Project name 10xTest must be snake_case and cannot start with a number or "pg_".')
  end


  it 'loads a project JSON definition' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/valid_project.json'
    )

    expect(project.valid?).to eq(true)
    expect(project.name).to eq('test')
    expect(project.name_full).to eq('Testing your immune system')

    expect(project.model_tree.map(&:name)).to eq([
      'project',
      'assay_pool',
      'document',
      'patient',
      'status',
      'symptom',
      'timepoint',
      'assay_name'
    ])
  end

  it 'reports errors for missing model information' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/missing_model_keys.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(8)
    expect(project.errors).to eq([
      "Missing required key for model assay_name: \"parent_link_type\".",
      "Missing required key for model document: \"parent_model_name\".",
      "Missing required key for model document: \"parent_link_type\".",
      "Parent model  for document does not exist in project.\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Missing required key for model assay_pool: \"identifier\".",
      "Missing required key for model patient: \"parent_model_name\".",
      "Parent model  for patient does not exist in project.\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Missing required key for model project: \"identifier\"."
    ])
  end

  it 'reports errors for invalid model information' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/invalid_models.json'
    )
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(11)
    expect(project.errors).to eq([
      "Model name assay_name!@ must be snake_case and can only consist of letters and \"_\".",
      "Invalid empty parent_model_name for model assay_name!@: \"\".",
      "Parent model  for assay_name!@ does not exist in project.\nCurrent models are [\"assay_name!@\", \"documents_2_keep\", \"assay_pool#\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Linked model, \"assay_pool\", on attribute assay_pool of model assay_name!@ does not exist!\nCurrent models are [\"assay_name!@\", \"documents_2_keep\", \"assay_pool#\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Model name documents_2_keep must be snake_case and can only consist of letters and \"_\".",
      "Invalid empty identifier for model documents_2_keep: \"\".",
      "Model name assay_pool# must be snake_case and can only consist of letters and \"_\".",
      "Invalid empty parent_link_type for model assay_pool#: \"\".",
      "Invalid parent_link_type for model assay_pool#: \"\".\nShould be one of [\"child\", \"collection\", \"table\"].",
      "Invalid parent_link_type for model patient: \"thingamajig\".\nShould be one of [\"child\", \"collection\", \"table\"].",
      "Invalid empty identifier for model project: \"\"."
    ])
  end

  it 'reports errors for invalid attribute information' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/attribute_errors.json'
    )

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(19)
    expect(project.errors).to eq([
      "Missing required key for model assay_name, attribute tube_name, validation: \"value\".",
      "Attribute name 123restricted in model assay_name must be snake_case and can only consist of letters, numbers, and \"_\".",
      "Missing required key for model assay_name, attribute assay_2_pool: \"desc\".",
      "Attribute name assay_2_pool in model assay_name should match the link_model_name, \"pool_deep_end\".",
      "Missing required key for model assay_name, attribute vendor: \"desc\".",
      "Invalid type for model assay_name, attribute vendor, validation: \"Lo que sea\".\nShould be one of [\"Array\", \"Range\", \"Regexp\"].",
      "Linked model, \"pool_deep_end\", on attribute assay_2_pool of model assay_name does not exist!\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Invalid empty attribute_type for model document, attribute document_desc: \"\".",
      "Invalid attribute_type for model document, attribute document_desc: \"\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Attribute name version!@ in model document must be snake_case and can only consist of letters, numbers, and \"_\".",
      "Invalid attribute_type for model document, attribute version!@: \"copacetic\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Invalid attribute_type for model document, attribute version_date: \"date\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Missing required key for model assay_pool, attribute biospecimen: \"attribute_type\".",
      "Missing required key for model assay_pool, attribute biospecimen: \"desc\".",
      "Invalid empty value for model assay_pool, attribute biospecimen, validation: \"\".",
      "Missing required key for model assay_pool, attribute cells_loaded: \"display_name\".",
      "Missing required key for model assay_pool, attribute project: \"desc\".",
      "Model \"assay_pool\" already belongs to parent model \"project\". Remove attribute \"project\".",
      "Missing required key for model timepoint, attribute day: \"desc\"."
    ])
  end

  it 'can get models as Magma::Models' do
    project = Etna::Clients::Magma::JsonProject.new(
      filepath: './spec/fixtures/create_project/valid_project.json'
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
    expect(assay_name.template.identifier).to eq('tube_name')
    expect(assay_name.template.parent).to eq('timepoint')

    status = source_models.model('status')
    expect(status.template.attributes.attribute('patient').attribute_type).to eq(Etna::Clients::Magma::AttributeType::PARENT)
    expect(status.template.attributes.attribute('patient').link_model_name).to eq('patient')
  end
end
