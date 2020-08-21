require 'webmock/rspec'
require 'json'
require_relative '../lib/commands'

class TestClass
  include Polyphemus::WithMetisWaiverHelpers

  def project
    :mvir1
  end
end

describe 'WithMetisWaiverHelpers Module' do
  let(:test_class) { TestClass.new }

  before(:each) do
    route_payload = JSON.generate([
      {:method=>"GET", :route=>"/:project_name/list_all_folders/:bucket_name", :name=>"folder_list_all_folders", :params=>["project_name", "bucket_name"]},
      {:method=>"GET", :route=>"/:project_name/list/:bucket_name/*folder_path", :name=>"folder_list", :params=>["project_name", "bucket_name", "folder_path"]},
      {:method=>"POST", :route=>"/:project_name/folder/rename/:bucket_name/*folder_path", :name=>"folder_rename", :params=>["project_name", "bucket_name", "folder_path"]},
      {:method=>"POST", :route=>"/:project_name/folder/create/:bucket_name/*folder_path", :name=>"folder_create", :params=>["project_name", "bucket_name", "folder_path"]}
    ])

    stub_request(:options, "https://metis.test/").
      to_return({
        status: 200,
        headers: {
          'Content-Type': 'application/json'
        },
        body: route_payload
      })
    stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
      .to_return({
        status: 200,
        headers: {
        'Content-Type' => 'application/json'
        },
        body: JSON.parse(File.read('spec/fixtures/metis_release_folder_fixture.json')).to_json
      })

    stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
      .to_return({
        status: 200,
        headers: {
        'Content-Type' => 'application/json'
        },
        body: JSON.parse(File.read('spec/fixtures/metis_restrict_folder_fixture.json')).to_json
      })
  end

  it 'fetches all data folders when called' do
    test_class.send('release_folders')
    test_class.send('restrict_folders')
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
  end

  it 'throws exception when does not find a valid patient to release' do
    stub_parent_exists

    expect {
      test_class.release_patient_data('Dan')
    }.to raise_error(Etna::Error)

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'attempts to rename all patient folders when releasing' do
    stub_parent_exists({status: 422, bucket: RELEASE_BUCKET})
    stub_create_folder({bucket: RELEASE_BUCKET})
    stub_rename_folder({bucket: RESTRICT_BUCKET})

    test_class.release_patient_data('Danielle')

    # There are two folders in assay/processed and two in assay/raw for each patient, in the fixtures
    # This should ignore the "summary" sub-folder that is under a patient
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}/).times(4)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/).times(4)
  end

  it 'throws exception when does not find a valid patient to restrict' do
    stub_parent_exists

    expect {
      test_class.restrict_patient_data('Danielle')
    }.to raise_error(Etna::Error)

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'attempts to rename all patient folders when restricting' do
    stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.restrict_patient_data('Dan')

    # There are two folders in assay/processed and two in assay/raw for each patient, in the fixtures
    # This should ignore the "summary" sub-folder that is under a patient
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/).times(4)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/).times(4)
  end

  it 'throws exception when does not find a valid pool to release' do
    expect {
      test_class.release_pool_data('pool-a')
    }.to raise_error(Etna::Error)

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/)
  end

  it 'attempts to rename all pool folders when releasing' do
    stub_parent_exists({status: 422, bucket: RELEASE_BUCKET})
    stub_create_folder({bucket: RELEASE_BUCKET})
    stub_rename_folder({bucket: RESTRICT_BUCKET})

    test_class.release_pool_data('pool-b')

    # There is one folder in assay/processed and one in assay/raw for each pool, in the fixtures
    # This should ignore the "summary" sub-folder that is under the pool
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}/).twice
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/).twice
  end

  it 'throws exception when does not find a valid pool to restrict' do
    expect {
      test_class.restrict_pool_data('pool-b')
    }.to raise_error(Etna::Error)

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'attempts to rename all pool folders when restricting' do
    stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.restrict_pool_data('pool-a')

    # There is one folder in assay/processed and one in assay/raw for each pool, in the fixtures
    # This should ignore the "summary" sub-folder that is under the pool
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/).twice
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/).twice
  end

  it 'throws exception if checking for parent results in non-422 error' do
    stub_parent_exists({status: 403, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    expect {
      test_class.restrict_pool_data('pool-a')
    }.to raise_error(Etna::Error)

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end
end