require 'webmock/rspec'
require 'json'
require_relative '../lib/etna/clients/metis'

describe 'Metis Client class' do
  let(:test_class) { Etna::Clients::Metis.new(token: 'fake-token', host: 'https://metis.test') }

  before(:each) do
    stub_metis_setup
  end

  it 'fetches all data folders when called' do
    test_class.list_all_folders(Etna::Clients::Metis::ListFoldersRequest.new(
      project_name: 'test',
      bucket_name: RELEASE_BUCKET
    ))
    test_class.list_all_folders(Etna::Clients::Metis::ListFoldersRequest.new(
      project_name: 'test',
      bucket_name: RESTRICT_BUCKET
    ))
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
  end

  it 'identifies the parent folder given a folder path' do
    expect(test_class.send('parent_folder_path', 'foo')).to eq('')
    expect(test_class.send('parent_folder_path', 'labors/monsters/avatars')).to eq('labors/monsters')
  end

  it 'creates a parent folder' do
    stub_create_folder({bucket: RELEASE_BUCKET})

    test_class.create_folder(Etna::Clients::Metis::CreateFolderRequest.new(
      project_name: 'test',
      bucket_name: RELEASE_BUCKET,
      folder_path: 'root_folder/sub_folder/grandchild_folder'
    ))

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'returns that root folders\' parent folder do exist' do
    stub_parent_exists({status: 200, bucket: RELEASE_BUCKET})

    expect(test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: ''))).to eq (true)

    expect(WebMock).not_to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}/)
  end

  it 'returns true if folder exists' do
    stub_parent_exists({status: 200, bucket: RELEASE_BUCKET})

    expect(test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/sub_folder'))).to eq (true)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'returns false if folder does not exist' do
    stub_parent_exists({status: 422, bucket: RELEASE_BUCKET})

    expect(test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/sub_folder'))).to eq(false)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder/)
  end

  it 'raises exception if folder exists check returns exception' do
    stub_parent_exists({status: 403, bucket: RELEASE_BUCKET})

    expect {
      test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
        bucket_name: RELEASE_BUCKET,
        project_name: 'test',
        folder_path: 'root_folder/sub_folder'))
    }.to raise_error(Etna::Error)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'renames a folder and creates parent if requested' do
    stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(Etna::Clients::Metis::RenameFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/subfolder',
      new_bucket_name: RESTRICT_BUCKET,
      new_folder_path: 'new_root_folder/new_subfolder',
      create_parent: true))

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a folder but does not create parent if not requested' do
    stub_parent_exists({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(Etna::Clients::Metis::RenameFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/subfolder',
      new_bucket_name: RESTRICT_BUCKET,
      new_folder_path: 'new_root_folder/new_subfolder'))

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a folder and does not create parent for new root folder' do
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(Etna::Clients::Metis::RenameFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder',
      new_bucket_name: RESTRICT_BUCKET,
      new_folder_path: 'new_root_folder',
      create_parent: true))

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/list\/#{RESTRICT_BUCKET}/)
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a folder and does not create parent when not needed' do
    stub_parent_exists({status: 200, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_folder({bucket: RELEASE_BUCKET})

    test_class.rename_folder(Etna::Clients::Metis::RenameFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/subfolder',
      new_bucket_name: RESTRICT_BUCKET,
      new_folder_path: 'new_root_folder/new_subfolder',
      create_parent: true))

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'fetches all data folders when called' do
    test_class.folders(project_name: PROJECT, bucket_name: RELEASE_BUCKET)
    test_class.folders(project_name: PROJECT, bucket_name: RESTRICT_BUCKET)
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RELEASE_BUCKET}/)
    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list_all_folders\/#{RESTRICT_BUCKET}/)
  end

  it 'can find folders and files' do
    stub_find({bucket: RELEASE_BUCKET})

    test_class.find(Etna::Clients::Metis::FindRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      params: [Etna::Clients::Metis::FindParam.new(
        attribute: 'name',
        predicate: 'glob',
        value: 'folder/fo*'
      )]))

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/find\/#{RELEASE_BUCKET}/).
      with(body: hash_including({
      params: [{
        attribute: 'name',
        predicate: 'glob',
        value: 'folder/fo*'
    }]}))
  end

  it 'can add multiple find params' do
    stub_find({bucket: RELEASE_BUCKET})

    find_request = Etna::Clients::Metis::FindRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      params: [])

    find_request.add_param(Etna::Clients::Metis::FindParam.new(
      attribute: 'name',
      predicate: 'glob',
      value: 'folder/fo*'
    ))

    test_class.find(find_request)

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/find\/#{RELEASE_BUCKET}/).
      with(body: hash_including({
      params: [{
        attribute: 'name',
        predicate: 'glob',
        value: 'folder/fo*'
    }]}))
  end
end