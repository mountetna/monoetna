require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::ProjectValidator do
  before(:each) do
  end

  def get_and_validate_project(filepath)
    project_json = JSON.parse(File.read(filepath))
    magma_project = Etna::Clients::Magma::Project.new(
      Etna::Clients::Magma::ConverterBase.convert_project_user_json_to_magma_json(project_json))

    project = Etna::Clients::Magma::ProjectValidator.new(magma_project)
    project.validate
    project
  end

  it 'reports errors for missing project keys' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/missing_project_keys.json')

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(3)
    expect(project.errors).to eq([
      'Missing required key for root project: "project_name".',
      'Missing required key for root project: "project_name_full".',
      'Project name  must be snake_case and cannot start with a number or "pg_".'
    ])
  end

  it 'reports errors for blank project keys' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/blank_project_keys.json')

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(3)
    expect(project.errors).to eq([
      'Invalid empty project_name for root project: "".',
      'Invalid empty project_name_full for root project: "".',
      'Project name  must be snake_case and cannot start with a number or "pg_".'
    ])
  end

  it 'reports errors for invalid project name' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/invalid_project_name.json')
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(1)
    expect(project.errors).to eq([
      'Project name 10xTest must be snake_case and cannot start with a number or "pg_".'])
  end

  it 'reports errors for missing model information' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/missing_model_keys.json')

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(8)
    expect(project.errors).to eq([
      "Parent model \"\" for document does not exist in project.\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Parent model \"\" for patient does not exist in project.\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Missing required key for model assay_name: \"parent_link_type\".",
      "Missing required key for model document: \"parent_model_name\".",
      "Missing required key for model document: \"parent_link_type\".",
      "Missing required key for model assay_pool: \"identifier\".",
      "Missing required key for model patient: \"parent_model_name\".",
      "Missing required key for model project: \"identifier\"."
    ])
  end

  it 'reports errors for invalid model information' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/invalid_models.json')
    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(14)
    expect(project.errors).to eq([
      "Parent model \"\" for assay_name!@ does not exist in project.\nCurrent models are [\"assay_name!@\", \"documents_2_keep\", \"assay_pool#\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Linked model, \"assay_pool\", on attribute assay_pool does not exist!\nCurrent models are [\"assay_name!@\", \"documents_2_keep\", \"assay_pool#\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Parent model \" project\" for patient does not exist in project.\nCurrent models are [\"assay_name!@\", \"documents_2_keep\", \"assay_pool#\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Model name assay_name!@ must be snake_case and can only consist of letters and \"_\".",
      "Invalid empty parent_model_name for model assay_name!@: \"\".",
      "Model name documents_2_keep must be snake_case and can only consist of letters and \"_\".",
      "Invalid empty identifier for model documents_2_keep: \"\".",
      "Model name assay_pool# must be snake_case and can only consist of letters and \"_\".",
      "Invalid empty parent_link_type for model assay_pool#: \"\".",
      "Invalid parent_link_type for model assay_pool#: \"\".\nShould be one of [\"child\", \"collection\", \"table\"].",
      "Invalid parent_link_type for model patient: \"thingamajig\".\nShould be one of [\"child\", \"collection\", \"table\"].",
      "Invalid empty identifier for model project: \"\".",
      "Invalid parent_link_type for model status: \"table \".\nShould be one of [\"child\", \"collection\", \"table\"].",
      "Missing required key for model status: \"identifier\"."
    ])
  end

  it 'reports errors for invalid attribute information' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/attribute_errors.json')

    expect(project.valid?).to eq(false)
    expect(project.errors.length).to eq(21)

    expect(project.errors).to eq([
      "Linked model, \"pool_deep_end\", on attribute assay_2_pool does not exist!\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Model \"assay_pool\" already belongs to parent model \"project\". Remove attribute \"project\".",
      "Missing required key for attribute tube_name, validation: \"value\".",
      "Attribute name \"123restricted\" must be snake_case and can only consist of letters, numbers, and \"_\".",
      "Missing required key for attribute assay_2_pool: \"desc\".",
      "Attribute name assay_2_pool should match the link_model_name, \"pool_deep_end\".",
      "Missing required key for attribute vendor: \"desc\".",
      "Invalid type for attribute vendor, validation: \"Lo que sea\".\nShould be one of [\"Array\", \"Range\", \"Regexp\"].",
      "Invalid empty attribute_type for attribute document_desc: \"\".",
      "Invalid attribute_type for attribute document_desc: \"\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Attribute name \"version!@\" must be snake_case and can only consist of letters, numbers, and \"_\".",
      "Invalid attribute_type for attribute version!@: \"copacetic\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Invalid attribute_type for attribute version_date: \"date\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Attribute key \"biospecimen \" must match attribute_name \"biospecimen\".",
      "Missing required key for attribute biospecimen: \"attribute_type\".",
      "Missing required key for attribute biospecimen: \"desc\".",
      "Invalid empty value for attribute biospecimen, validation: \"\".",
      "Attribute key \"cells_loaded \" must match attribute_name \"cells_loaded\".",
      "Missing required key for attribute cells_loaded: \"display_name\".",
      "Missing required key for attribute project: \"desc\".",
      "Missing required key for attribute day: \"desc\"."
    ])
  end

  it 'reports reports no errors for valid project JSON' do
    project = get_and_validate_project(
      './spec/fixtures/create_project/valid_project.json')

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
    expect(validator.errors.length).to eq(7)

    expect(validator.errors).to eq([
      "Missing required key for attribute notes: \"attribute_type\".",
      "Missing required key for attribute notes: \"display_name\".",
      "Missing required key for attribute notes: \"desc\".",
      "Attribute \"notes\" does not exist in model assay_name.",
      "Attribute \"notes\" does not exist in model assay_name.",
      "Attribute \"assay_pool\" already exists in model assay_name.",
      "Attribute \"assay_name\" already exists in model assay_pool."
    ])
  end

  it 'reports errors for invalid add_attribute actions' do
    validator = get_and_validate_actions(
      './spec/fixtures/attribute_actions/add_attribute_invalid.json')

    expect(validator.valid?).to eq(false)
    expect(validator.errors.length).to eq(4)
    expect(validator.errors).to eq([
      "Missing required key for attribute notes: \"attribute_type\".",
      "Missing required key for attribute notes: \"desc\".",
      "Model \"asay_name\" does not exist in project.",
      "Attribute \"vendor\" already exists in model assay_name."])
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
