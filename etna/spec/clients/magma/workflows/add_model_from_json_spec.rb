require 'webmock/rspec'
require 'json'

describe Etna::Clients::Magma::AddModelFromJsonWorkflow do
  let(:magma_client) {Etna::Clients::Magma.new(
    token: '123',
    host: MAGMA_HOST)}

  before(:each) do
    @all_updates = []
    stub_magma_models({
      models: {
      assay_name: model_stamp('assay_name'),
      assay_pool: model_stamp('assay_pool'),
      project: model_stamp('project'),
      timepoint: model_stamp('timepoint'),
      patient: model_stamp('patient'),
      document: model_stamp('document'),
      status: model_stamp('status'),
      symptom: model_stamp('symptom')
    }})
    stub_magma_update_model
    stub_magma_update
  end

  def updates
    @all_updates.inject({}) do |acc, n|
      n.keys.each do |k|
        (acc[k] ||= {}).update(n[k])
      end
      acc
    end
  end

  it 'adds the model to magma' do
    # Need different responses for /retrieve for ensure_model_tree
    stub_request(:post, /#{MAGMA_HOST}\/retrieve/)
      .to_return({body: {
        models: {
        assay_name: model_stamp('assay_name'),
        assay_pool: model_stamp('assay_pool'),
        project: model_stamp('project'),
        timepoint: model_stamp('timepoint'),
        patient: model_stamp('patient'),
        document: model_stamp('document'),
        status: model_stamp('status'),
        symptom: model_stamp('symptom')
      }}.to_json}).then
      .to_return({body: {
        models: {
        assay_name: model_stamp('assay_name'),
        assay_pool: model_stamp('assay_pool'),
        project: model_stamp('project'),
        timepoint: model_stamp('timepoint'),
        patient: model_stamp('patient'),
        document: model_stamp('document'),
        status: model_stamp('status'),
        symptom: model_stamp('symptom'),
        assay_two: model_stamp('assay_two')
      }}.to_json})

    workflow = Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
      magma_client: magma_client,
      project_name: PROJECT,
      model_name: 'assay_two',
      filepath: './spec/fixtures/add_model/add_model_fixture_valid.json'
    )
    workflow.add!

    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna 123'}).
      with { |req| req.body.include?('add_model') }.times(1)
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna 123'}).
      with { |req| req.body.include?('add_attribute') }.times(3)  # 3 attributes in the fixture
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna 123'}).
      with { |req| req.body.include?('add_link') }.times(2)
    # This is for the link ... this gets called twice because our Webmock
    #   does not get updated after the first call. Should not happen in a
    #   real-world environment.


    # Make sure the assay_two identifier validation is submitted.
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update_model/).
      with(headers: {Authorization: 'Etna 123'}).
      with { |req| req.body.include?('update_attribute') }
  end

  it 'raises an exception if the model already exists' do
    expected_msg = %{Model JSON has errors:
  * Model assay_name already exists in project #{PROJECT}!}

    expect {
      Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
        magma_client: magma_client,
        project_name: PROJECT,
        model_name: 'assay_name',
        filepath: './spec/fixtures/add_model/add_model_fixture_valid.json'
      )
    }.to raise_error(Exception, expected_msg)
  end

  it 'raises exception for an invalid model JSON' do

    expected_msg = %{Model JSON has errors:
  * Model name assay2 must be snake_case and can only consist of letters and "_".
  * Invalid parent_link_type for model assay2: \"magma\".
\tShould be one of [\"child\", \"collection\", \"table\"].
  * Invalid type for model assay2, attribute vendor, validation: \"array\".
\tShould be one of [\"Array\", \"Range\", \"Regexp\"].
  * Parent model paper_airplanes for assay2 does not exist in project.
\tCurrent models are [\"assay_name\", \"assay_pool\", \"project\", \"timepoint\", \"patient\", \"document\", \"status\", \"symptom\"].
  * Linked model, "assay2_pool", on attribute assay2_pool of model assay2 does not exist!
\tCurrent models are [\"assay_name\", \"assay_pool\", \"project\", \"timepoint\", \"patient\", \"document\", \"status\", \"symptom\"].}

    expect {
      Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
        magma_client: magma_client,
        project_name: PROJECT,
        model_name: 'assay2',
        filepath: './spec/fixtures/add_model/add_model_fixture_invalid.json'
      )
    }.to raise_error(Exception, expected_msg)
  end
end
