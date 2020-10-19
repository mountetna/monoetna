require 'webmock/rspec'
require 'json'
require 'pry'

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
      magma_json.dig('models', 'assay_name', 'template', 'name')).to eq('assay_name')
  end

  it 'correctly assigns the identifier type' do
    project_json = JSON.parse(File.read('./spec/fixtures/create_project/valid_project.json'))

    magma_json = Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json)

    expect(
      magma_json.dig('models', 'assay_name', 'template', 'attributes', 'tube_name', 'type')).to eq('identifier')
  end
end
