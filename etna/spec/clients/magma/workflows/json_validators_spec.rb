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
    expect(project.errors.first).to eq('Project name 10xTest must be snake_case and cannot start with a number or "pg_".')
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
    expect(project.errors.length).to eq(19)
    expect(project.errors).to eq([
      "Linked model, \"pool_deep_end\", on attribute assay_2_pool does not exist!\nCurrent models are [\"assay_name\", \"document\", \"assay_pool\", \"patient\", \"project\", \"status\", \"symptom\", \"timepoint\"].",
      "Model \"assay_pool\" already belongs to parent model \"project\". Remove attribute \"project\".",
      "Missing required key for attribute tube_name, validation: \"value\".",
      "Attribute name 123restricted must be snake_case and can only consist of letters, numbers, and \"_\".",
      "Missing required key for attribute assay_2_pool: \"desc\".",
      "Attribute name assay_2_pool should match the link_model_name, \"pool_deep_end\".",
      "Missing required key for attribute vendor: \"desc\".",
      "Invalid type for attribute vendor, validation: \"Lo que sea\".\nShould be one of [\"Array\", \"Range\", \"Regexp\"].",
      "Invalid empty attribute_type for attribute document_desc: \"\".",
      "Invalid attribute_type for attribute document_desc: \"\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Attribute name version!@ must be snake_case and can only consist of letters, numbers, and \"_\".",
      "Invalid attribute_type for attribute version!@: \"copacetic\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Invalid attribute_type for attribute version_date: \"date\".\nShould be one of [\"boolean\", \"date_time\", \"file\", \"float\", \"image\", \"integer\", \"link\", \"match\", \"matrix\", \"string\", \"table\"].",
      "Missing required key for attribute biospecimen: \"attribute_type\".",
      "Missing required key for attribute biospecimen: \"desc\".",
      "Invalid empty value for attribute biospecimen, validation: \"\".",
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
