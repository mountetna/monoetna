require "webmock/rspec"
require "json"
require_relative "../lib/commands"
require_relative "../lib/metis/mvir1_waiver"

describe "Mvir1Waiver class" do
  let(:test_class) {
    Mvir1Waiver.new(
      metis_client: Etna::Clients::Metis.new(token: "fake-token", host: "https://metis.test"),
      project_name: PROJECT,
    )
  }

  def stub_no_data
    stub_bucket_find(
      bucket: RELEASE_BUCKET,
      response_body: {
        folders: [],
      },
    )
    stub_bucket_find(
      bucket: RESTRICT_BUCKET,
      response_body: {
        folders: [],
      },
    )
  end

  before(:each) do
    stub_metis_setup
  end

  it "does not throw an exception when does not find a valid patient to release" do
    stub_parent_exists
    stub_no_data
    test_class.release_patient_data("Dan")
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it "attempts to rename all patient folders when releasing" do
    stub_parent_exists({ status: 422, bucket: RELEASE_BUCKET })
    stub_create_folder({ bucket: RELEASE_BUCKET })
    stub_rename_folder({ bucket: RESTRICT_BUCKET })
    stub_bucket_find(
      bucket: RESTRICT_BUCKET,
      response_body: {
        folders: [
          create_metis_folder("Danielle-D7-ASSAY1", "assay/processed/Danielle-D7-ASSAY1").raw,
          create_metis_folder("Danielle-D14-ASSAY1", "assay/processed/Danielle-D14-ASSAY1").raw,
          create_metis_folder("Danielle-D7", "assay/raw/Danielle-D7").raw,
          create_metis_folder("Danielle-D14-ASSAY2", "assay/raw/Danielle-D14-ASSAY2").raw,
        ],
      },
    )
    stub_bucket_find(
      bucket: RELEASE_BUCKET,
      response_body: {
        folders: [],
      },
    )

    test_class.release_patient_data("Danielle")

    # There are two folders in assay/processed and two in assay/raw for each patient, in the fixtures
    # This should ignore the "summary" sub-folder that is under a patient
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}/).times(4)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/).times(4)
  end

  it "does not throw  an exception when does not find a valid patient to restrict" do
    stub_parent_exists
    stub_no_data
    test_class.restrict_patient_data("Danielle")
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it "attempts to rename all patient folders when restricting" do
    stub_parent_exists({ status: 422, bucket: RESTRICT_BUCKET })
    stub_create_folder({ bucket: RESTRICT_BUCKET })
    stub_rename_folder({ bucket: RELEASE_BUCKET })
    stub_bucket_find(
      bucket: RELEASE_BUCKET,
      response_body: {
        folders: [
          create_metis_folder("Dan-D0-ASSAY1", "assay/processed/Dan-D0-ASSAY1").raw,
          create_metis_folder("Dan-D12-ASSAY1", "assay/processed/Dan-D12-ASSAY1").raw,
          create_metis_folder("Dan-D0-ASSAY3", "assay/raw/Dan-D0-ASSAY3").raw,
          create_metis_folder("Dan-D12-ASSAY3", "assay/raw/Dan-D12-ASSAY3").raw,
        ],
      },
    )
    stub_bucket_find(
      bucket: RESTRICT_BUCKET,
      response_body: {
        folders: [],
      },
    )

    test_class.restrict_patient_data("Dan")

    # There are two folders in assay/processed and two in assay/raw for each patient, in the fixtures
    # This should ignore the "summary" sub-folder that is under a patient
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/).times(4)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/).times(4)
  end

  it "does not throw an exception when does not find a valid pool to release" do
    stub_no_data
    test_class.release_pool_data("pool-a")
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/)
  end

  it "attempts to rename all pool folders when releasing" do
    stub_parent_exists({ status: 422, bucket: RELEASE_BUCKET })
    stub_create_folder({ bucket: RELEASE_BUCKET })
    stub_rename_folder({ bucket: RESTRICT_BUCKET })
    stub_bucket_find(
      bucket: RESTRICT_BUCKET,
      response_body: {
        folders: [
          create_metis_folder("pool-b", "assay/processed/pool-b").raw,
          create_metis_folder("pool-b", "assay/raw/pool-b").raw,
        ],
      },
    )
    stub_bucket_find(
      bucket: RELEASE_BUCKET,
      response_body: {
        folders: [],
      },
    )

    test_class.release_pool_data("pool-b")

    # There is one folder in assay/processed and one in assay/raw for each pool, in the fixtures
    # This should ignore the "summary" sub-folder that is under the pool
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}/).twice
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RESTRICT_BUCKET}/).twice
  end

  it "does not throw an exception when does not find a valid pool to restrict" do
    stub_no_data
    test_class.restrict_pool_data("pool-b")
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it "attempts to rename all pool folders when restricting" do
    stub_parent_exists({ status: 422, bucket: RESTRICT_BUCKET })
    stub_create_folder({ bucket: RESTRICT_BUCKET })
    stub_rename_folder({ bucket: RELEASE_BUCKET })
    stub_bucket_find(
      bucket: RELEASE_BUCKET,
      response_body: {
        folders: [
          create_metis_folder("pool-a", "assay/processed/pool-a").raw,
          create_metis_folder("pool-a", "assay/raw/pool-a").raw,
        ],
      },
    )
    stub_bucket_find(
      bucket: RESTRICT_BUCKET,
      response_body: {
        folders: [],
      },
    )

    test_class.restrict_pool_data("pool-a")

    # There is one folder in assay/processed and one in assay/raw for each pool, in the fixtures
    # This should ignore the "summary" sub-folder that is under the pool
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/).twice
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/).twice
  end

  it "throws exception if checking for parent results in non-422 error" do
    stub_parent_exists({ status: 403, bucket: RESTRICT_BUCKET })
    stub_create_folder({ bucket: RESTRICT_BUCKET })
    stub_rename_folder({ bucket: RELEASE_BUCKET })
    stub_bucket_find(
      bucket: RELEASE_BUCKET,
      response_body: {
        folders: [
          create_metis_folder("pool-a", "assay/processed/pool-a").raw,
          create_metis_folder("pool-a", "assay/raw/pool-a").raw,
        ],
      },
    )
    stub_bucket_find(
      bucket: RESTRICT_BUCKET,
      response_body: {
        folders: [],
      },
    )

    expect {
      test_class.restrict_pool_data("pool-a")
    }.to raise_error(Etna::Error)

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end
end
