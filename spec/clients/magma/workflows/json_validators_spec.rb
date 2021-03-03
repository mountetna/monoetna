require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::ProjectValidator do
  before(:each) do
  end

  def get_and_validate_project(**args)
    project = Etna::Clients::Magma::ProjectValidator.new(**args)
    project.validate
    project
  end

  it 'reports errors for missing project keys' do
    project = get_and_validate_project

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(3)
    expect(project.errors).to eq([
      'Missing required key for root project: "project_name".',
      'Missing required key for root project: "project_name_full".',
      'Project name  must be snake_case and cannot start with a number or "pg_".'
    ])
  end

  it 'reports errors for blank project keys' do
    project = get_and_validate_project(project_name: '', project_name_full: '')
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(3)
    expect(project.errors).to eq([
      'Invalid empty project_name for root project: "".',
      'Invalid empty project_name_full for root project: "".',
      'Project name  must be snake_case and cannot start with a number or "pg_".'
    ])
  end

  it 'reports errors for invalid project name' do
    project = get_and_validate_project(project_name: '10xTest', project_name_full: 'abc')
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(1)
    expect(project.errors).to eq([
      'Project name 10xTest must be snake_case and cannot start with a number or "pg_".'])
  end

  it 'reports reports no errors for valid project JSON' do
    project = get_and_validate_project(project_name: 'myproject', project_name_full: 'My Project')

    expect(project.valid?).to eq(true)
  end
end


describe Etna::Clients::Magma::AttributeActionsValidator do
  let(:project_models) {
    project_json = JSON.parse(File.read(
      './spec/fixtures/attribute_actions/test_project_magma_models.json'
    ))
    Etna::Clients::Magma::Project.new(project_json).models
  }

  before(:each) do
  end

  def get_and_validate_actions(filepath)
    actions_json = JSON.parse(File.read(filepath))

    converter = Etna::Clients::Magma::AttributeActionsConverter.new(actions_json)
    actions = converter.convert

    validator = Etna::Clients::Magma::AttributeActionsValidator.new(
      actions,
      project_models)
    validator.validate
    validator
  end

  it 'reports no errors with multiple valid actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/multiple_actions_valid.json')

    expect(validator.valid?).to eq(true)
  end

  it 'reports errors with multiple invalid actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/multiple_actions_invalid.json')

    expect(validator.valid?).to eq(false)
    expect(validator.errors.length).to eq(5)
  end

  it 'reports errors for invalid add_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/add_attribute_invalid.json')

    expect(validator.valid?).to eq(false)
    expect(validator.errors.length).to eq(3)
  end

  it 'reports no errors for valid add_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/add_attribute_valid.json')

    expect(validator.valid?).to eq(true)
  end

  it 'reports errors for invalid add_link actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/add_link_invalid.json')

    expect(validator.valid?).to eq(false)
    expect(validator.errors.length).to eq(8)
    expect(validator.errors).to eq([
      "Links {:model_name=>\"assay_name\", :attribute_name=>\"document\", :type=>\"link\"} and {:model_name=>\"patient\", :attribute_name=>\"assay_name\", :type=>\"collection\"} must point to each other.",
      "Model \"patient_demographics\" does not exist in project.",
      "Links {:model_name=>\"assay_name\", :attribute_name=>\"document_two\", :type=>\"link\"} and {:model_name=>\"patient_demographics\", :attribute_name=>\"assay_name_too\", :type=>\"collection\"} must point to each other.",
      "You must have one \"link\" and one \"collection\" type in the links.",
      "Must include two link entries, each with \"model_name\", \"attribute_name\", and \"type\".",
      "Missing required key for link {:model_name=>\"assay_name\", :attribute_name=>\"document\", :type=>nil}: \"type\".",
      "Missing required key for link {:model_name=>\"assay_name\", :attribute_name=>\"document\", :type=>nil}: \"type\".",
      "You must have one \"link\" and one \"collection\" type in the links."
    ])
  end

  it 'reports no errors for valid add_link actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/add_link_valid.json')

    expect(validator.valid?).to eq(true)
  end

  it 'reports errors for invalid rename_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/rename_attribute_invalid.json')

    expect(validator.valid?).to eq(false)
    expect(validator.errors.length).to eq(2)
    expect(validator.errors).to eq([
      "Attribute \"notes\" does not exist in model assay_name.",
      "Model \"asay_name\" does not exist in project."
    ])
  end

  it 'reports errors for invalid rename_attribute keyword' do
    expected_message = %{Exception while parsing {:action_name=>"rename_attribute", :model_name=>"assay_name", :old_attribute_name=>"vendor", :attribute_name=>"vendor_2"}.
unknown keywords: old_attribute_name}

    expect {
      get_and_validate_actions(
        './spec/fixtures/attribute_actions/rename_attribute_invalid_keyword.json')
    }.to raise_error(ArgumentError, expected_message)
  end

  it 'reports no errors for valid rename_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/rename_attribute_valid.json')

    expect(validator.valid?).to eq(true)
  end

  it 'reports errors for invalid update_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/update_attribute_invalid.json')

    expect(validator.valid?).to eq(false)
    expect(validator.errors.length).to eq(2)
    expect(validator.errors).to eq([
      "Attribute \"notes\" does not exist in model assay_name.",
      "Model \"asay_name\" does not exist in project."
    ])
  end

  it 'reports no errors for valid update_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/update_attribute_valid.json')

    expect(validator.valid?).to eq(true)
  end
end
