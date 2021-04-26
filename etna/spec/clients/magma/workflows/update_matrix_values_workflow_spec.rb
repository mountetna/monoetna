require "webmock/rspec"
require "json"

describe Etna::Clients::Magma::UpdateMatrixValuesWorkflow do
  let(:magma_client) {
    Etna::Clients::Magma.new(
      token: TEST_TOKEN,
      host: MAGMA_HOST,
    )
  }

  before(:each) do
    stub_magma_models(
      JSON.parse(File.read("./spec/fixtures/magma/magma_test_model.json"))
    )
    stub_magma_update_json
  end

  it "inserts zero for values that are not provided" do
    model_name = "model_two"
    attribute_name = "contributions"
    record_name = "record_1_matrix"

    workflow = Etna::Clients::Magma::UpdateMatrixValuesWorkflow.new(
      magma_client: magma_client,
      project_name: PROJECT,
      model_name: model_name,
      attribute_name: attribute_name,
      filepath: "./spec/fixtures/magma/#{record_name}_counts.csv",
      execute: true,
    )

    workflow.upload_values

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                         .with(body: hash_including({
                                 "revisions": {
                                   "#{model_name}": {
                                     "#{record_name}": {
                                       "#{attribute_name}": [3, 2, 0, 92],
                                     },
                                   },
                                 },
                                 "project_name": PROJECT,
                               }))
  end

  it "inserts zero for values that are not provided with TSV" do
    model_name = "model_two"
    attribute_name = "contributions"
    record_name = "record_1_matrix"

    workflow = Etna::Clients::Magma::UpdateMatrixValuesWorkflow.new(
      magma_client: magma_client,
      project_name: PROJECT,
      model_name: model_name,
      attribute_name: attribute_name,
      filepath: "./spec/fixtures/magma/#{record_name}_counts.tsv",
      execute: true,
      tsv: true
    )

    workflow.upload_values

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                         .with(body: hash_including({
                                 "revisions": {
                                   "#{model_name}": {
                                     "#{record_name}": {
                                       "#{attribute_name}": [3, 2, 0, 92],
                                     },
                                   },
                                 },
                                 "project_name": PROJECT,
                               }))
  end
end
