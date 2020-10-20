require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::ConverterBase do
  before(:each) do
  end

  it 'preserves the project_name and project_name_full keys' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(magma_json['project_name']).to eq('test')
    expect(magma_json['project_name_full']).to eq('Testing your immune system')
  end

  it 'correctly injects the template key into the JSON for models' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    expect(project_json['models']['assay_name'].keys.include?('template')).to eq(false)

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(magma_json['models']['assay_name'].keys.include?('template')).to eq(true)
  end

  it 'correctly moves the model attributes to under the template key' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    expect(project_json['models']['assay_name'].keys.include?('attributes')).to eq(true)

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(magma_json['models']['assay_name']['template'].keys.include?('attributes')).to eq(true)
  end

  it 'correctly copies the model name to the template' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(
      magma_json.dig('models', 'assay_name', 'template', 'name')
    ).to eq('assay_name')
  end

  it 'correctly copies the model identifier to the template' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(
      magma_json.dig('models', 'assay_name', 'template', 'identifier')
    ).to eq('tube_name')
  end

  it 'correctly copies the model parent to the template' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(
      magma_json.dig('models', 'assay_name', 'template', 'parent')
    ).to eq('timepoint')
  end

  it 'correctly assigns the identifier type' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(
      magma_json.dig('models', 'assay_name', 'template', 'attributes', 'tube_name', 'attribute_type')
    ).to eq(Etna::Clients::Magma::AttributeType::IDENTIFIER)
  end
end

describe Etna::Clients::Magma::ProjectConverter do
  before(:each) do
  end

  def get_project(filepath)
    project_json = JSON.parse(File.read(filepath))
    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)
    Etna::Clients::Magma::Project.new(magma_json)
  end

  it 'loads a project JSON definition' do
    project = get_project(
      './spec/fixtures/create_project/valid_project.json')
    converter = Etna::Clients::Magma::ProjectConverter.new(project)

    expect(converter.model_tree.map(&:name)).to eq([
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

  it 'correctly adds missing attributes for Magma format' do
    proj = get_project(
      './spec/fixtures/create_project/valid_project.json')

    converter = Etna::Clients::Magma::ProjectConverter.new(proj)
    converter.convert!

    source_models = proj.models

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

describe Etna::Clients::Magma::AttributeActionsConverter do
  before(:each) do
  end

  it 'converts JSON attribute actions to Magma models' do
    actions_json = JSON.parse(File.read('./spec/fixtures/attribute_actions/multiple_actions_valid.json'))
    converter = Etna::Clients::Magma::AttributeActionsConverter.new(actions_json)
    actions = converter.convert
    expect(actions.map(&:class)).to eq([
      Etna::Clients::Magma::AddAttributeAction,
      Etna::Clients::Magma::UpdateAttributeAction,
      Etna::Clients::Magma::RenameAttributeAction,
      Etna::Clients::Magma::AddLinkAction
    ])
    expect(actions.first.type).to eq(actions_json.first['attribute_type'])
    expect(actions.last.links.first[:type]).to eq(actions_json.last['links'].first['attribute_type'])
  end
end