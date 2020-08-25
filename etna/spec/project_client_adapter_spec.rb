require 'webmock/rspec'
require 'json'
require_relative '../lib/etna/clients/metis'


describe 'Metis ProjectClientAdapter class' do
  let(:test_class) { Etna::Clients::ProjectClientAdapter.new(
      PROJECT,
      Etna::Clients::Metis.new(token: 'fake-token', host: 'https://metis.test')
  ) }

  before(:each) do
    stub_metis_setup
  end

  it 'fetches all data folders when called' do
    test_class.fetch_folders(RELEASE_BUCKET)
    test_class.fetch_folders(RESTRICT_BUCKET)
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
  end

  it 'identifies the parent folder given a folder path' do
    expect(test_class.parent_folder_path('foo')).to eq('')
    expect(test_class.parent_folder_path('labors/monsters/avatars')).to eq('labors/monsters')
  end

  it 'creates a parent folder' do
    stub_create_folder({bucket: RELEASE_BUCKET})

    test_class.create_parent_folder(RELEASE_BUCKET, 'root_folder/sub_folder/grandchild_folder')

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'returns that root folders\' parent folder do exist' do
    stub_parent_exists({status: 200, bucket: RELEASE_BUCKET})

    expect(test_class.parent_exists?(RELEASE_BUCKET, 'root_folder')).to eq (true)

    expect(WebMock).not_to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}/)
  end

  it 'returns true if subdirectories do have a parent folder' do
    stub_parent_exists({status: 200, bucket: RELEASE_BUCKET})

    expect(test_class.parent_exists?(
      RELEASE_BUCKET, 'root_folder/sub_folder/grandchild_folder')).to eq (true)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'returns false if subdirectories do not have a parent folder' do
    stub_parent_exists({status: 422, bucket: RELEASE_BUCKET})

    expect(test_class.parent_exists?(
      RELEASE_BUCKET, 'root_folder/sub_folder/grandchild_folder')).to eq(false)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder/)
  end

  it 'raises exception if parent folder check returns exception' do
    stub_parent_exists({status: 403, bucket: RELEASE_BUCKET})

    expect {
      test_class.parent_exists?(RELEASE_BUCKET, 'root_folder/sub_folder/grandchild_folder')
    }.to raise_error(Etna::Error)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'renames a folder and creates parent if needed' do
    stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(RELEASE_BUCKET, 'root_folder/subfolder', RESTRICT_BUCKET, 'new_root_folder/new_subfolder')

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a folder and does not create parent for new root folder' do
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(RELEASE_BUCKET, 'root_folder', RESTRICT_BUCKET, 'new_root_folder')

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/list\/#{RESTRICT_BUCKET}/)
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a folder and does not create parent when not needed' do
    stub_parent_exists({status: 200, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(RELEASE_BUCKET, 'root_folder/sub_folder', RESTRICT_BUCKET, 'new_root_folder/new_sub_folder')

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end
end