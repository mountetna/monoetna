require 'webmock/rspec'
require 'json'
require 'pry'

describe Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow do
  let(:magma_client) {Etna::Clients::Magma.new(
      token: '123',
      host: MAGMA_HOST)}
  let(:magma_crud) {Etna::Clients::Magma::MagmaCrudWorkflow.new(
      magma_client: magma_client, project_name: PROJECT)}

  before(:each) do
    stub_magma_models(
      JSON.parse(File.read('./spec/fixtures/magma/magma_test_model.json')))
    stub_magma_update
  end

  it 'raises exception for invalid models' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/magma_update_attributes_invalid_model.csv'
    )

    expect {
      workflow.update_attributes
    }.to raise_error(RuntimeError, 'Invalid model fake_model for project test.')

    expect(WebMock).not_to have_requested(:post, /#{MAGMA_HOST}\/update/)
  end

  it 'sends valid revisions to magma' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/magma_update_attributes_valid.csv'
    )

    workflow.update_attributes

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_two": {"234": {"height": 64.2,"invisible": true}}}
      }))
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_two": {"123": {"name": "Record #123","strength": 2,"invisible": true}}}
      }))
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_one": {"xyz": {"name": "First model of the day"}}}
      }))
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_one": {"098": {"name": nil}}}
      }))
  end

  it 'raises exception for invalid attribute name' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/magma_update_attributes_invalid_attribute.csv'
    )

    expect {
      workflow.update_attributes
    }.to raise_error(RuntimeError, 'Invalid attribute weight for model model_two.')

    expect(WebMock).not_to have_requested(:post, /#{MAGMA_HOST}\/update/)
  end

  it 'treats malformed boolean values as false' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/magma_update_attributes_invalid_boolean.csv'
    )

    workflow.update_attributes

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_two": {"234": {"invisible": false}}}
      }))
  end

  it 'casts values for integer attributes' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/magma_update_attributes_invalid_integer.csv'
    )

    workflow.update_attributes

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_two": {"234": {"strength": 3}}}
      }))
  end

  it 'casts values for float values' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/magma_update_attributes_invalid_float.csv'
    )

    workflow.update_attributes

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).
      with(body: hash_including({
        "revisions": {"model_two": {"234": {"strength": 0.0}}}
      }))
  end

  it 'raises exception for invalid filename' do
    workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflow.new(
      magma_crud: magma_crud,
      project_name: PROJECT,
      input_filename: './spec/fixtures/magma/nonexistent_input.csv'
    )

    expect {
      workflow.update_attributes
    }.to raise_error(StandardError)

    expect(WebMock).not_to have_requested(:post, /#{MAGMA_HOST}\/update/)
  end
end