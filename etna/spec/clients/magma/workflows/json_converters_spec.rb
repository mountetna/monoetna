require 'webmock/rspec'
require 'json'

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