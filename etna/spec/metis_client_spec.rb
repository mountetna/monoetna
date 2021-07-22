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
    stub_list_bucket({status: 200, bucket: RELEASE_BUCKET})

    expect(test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: ''))).to eq (true)

    expect(WebMock).not_to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}/)
  end

  it 'returns true if folder exists' do
    stub_list_bucket({status: 200, bucket: RELEASE_BUCKET})

    expect(test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/sub_folder'))).to eq (true)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'returns false if folder does not exist' do
    stub_list_bucket({status: 422, bucket: RELEASE_BUCKET})

    expect(test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      folder_path: 'root_folder/sub_folder'))).to eq(false)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder/)
  end

  it 'raises exception if folder exists check returns exception' do
    stub_list_bucket({status: 403, bucket: RELEASE_BUCKET})

    expect {
      test_class.folder_exists?(Etna::Clients::Metis::CreateFolderRequest.new(
        bucket_name: RELEASE_BUCKET,
        project_name: 'test',
        folder_path: 'root_folder/sub_folder'))
    }.to raise_error(Etna::Error)

    expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\/root_folder\/sub_folder?/)
  end

  it 'renames a folder and creates parent if requested' do
    stub_list_bucket({status: 422, bucket: RESTRICT_BUCKET})
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
    stub_list_bucket({status: 422, bucket: RESTRICT_BUCKET})
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
    stub_list_bucket({status: 200, bucket: RESTRICT_BUCKET})
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

  it 'renames a file and creates parent if requested' do
    stub_list_bucket({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_file({bucket: RELEASE_BUCKET})

    test_class.rename_file(Etna::Clients::Metis::RenameFileRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      file_path: 'root_folder/subfolder/myfile.txt',
      new_bucket_name: RESTRICT_BUCKET,
      new_file_path: 'new_root_folder/new_subfolder/yourfile.txt',
      create_parent: true))

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/file\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a file but does not create parent if not requested' do
    stub_list_bucket({status: 422, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_file({bucket: RELEASE_BUCKET})

    test_class.rename_file(Etna::Clients::Metis::RenameFileRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      file_path: 'root_folder/subfolder/myfile.txt',
      new_bucket_name: RESTRICT_BUCKET,
      new_file_path: 'new_root_folder/new_subfolder/yourfile.txt'))

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/file\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a file and does not create parent for root folder' do
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_file({bucket: RELEASE_BUCKET})

    test_class.rename_file(Etna::Clients::Metis::RenameFileRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      file_path: 'root_folder/myfile.txt',
      new_bucket_name: RESTRICT_BUCKET,
      new_file_path: 'yourfile.txt',
      create_parent: true))

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/list\/#{RESTRICT_BUCKET}/)
    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/file\/rename\/#{RELEASE_BUCKET}/)
  end

  it 'renames a file and does not create parent when not needed' do
    stub_list_bucket({status: 200, bucket: RESTRICT_BUCKET})
    stub_create_folder({bucket: RESTRICT_BUCKET})
    stub_rename_file({bucket: RELEASE_BUCKET})

    test_class.rename_file(Etna::Clients::Metis::RenameFileRequest.new(
      bucket_name: RELEASE_BUCKET,
      project_name: 'test',
      file_path: 'root_folder/myfile.txt',
      new_bucket_name: RESTRICT_BUCKET,
      new_file_path: 'yourfile.txt',
      create_parent: true))

    expect(WebMock).not_to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)
    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/file\/rename\/#{RELEASE_BUCKET}/)
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

  it 'can copy files' do
    stub_copy

    copy_request = Etna::Clients::Metis::CopyFilesRequest.new(
      project_name: 'test',
      revisions: [])

    copy_request.add_revision(Etna::Clients::Metis::CopyRevision.new(
      source: 'metis://foo',
      dest: 'metis://bar'
    ))

    test_class.copy_files(copy_request)

    expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/files\/copy/).
      with(body: hash_including({
      revisions: [{
        source: 'metis://foo',
        dest: 'metis://bar'
    }]}))
  end

  describe 'renames folders by regex' do
    it 'when dest folder does not exist' do
      stub_list_bucket({status: 422, bucket: RESTRICT_BUCKET})
      stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\//)
        .to_return({
          status: 200,
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [{
              "folder_name": "premessa",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "folder_path": "assay/processed/sample-10/premessa"
            }],
            files: [{
              "file_name": "data.txt",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "file_path": "assay/processed/sample-10/data.txt"
            }]
          }.to_json})
      stub_create_folder({bucket: RESTRICT_BUCKET})
      stub_rename_folder({bucket: RELEASE_BUCKET})
      stub_delete_folder({bucket: RELEASE_BUCKET})

      source_folders = Etna::Clients::Metis::Folders.new([{
        "folder_name": "sample-42",
        "bucket_name": RELEASE_BUCKET,
        "project_name": PROJECT,
        "folder_path": "assay/processed/sample-42"
      }])

      test_class.rename_folders_by_regex(
        source_bucket: RELEASE_BUCKET,
        project_name: PROJECT,
        dest_bucket: RESTRICT_BUCKET,
        source_folders: source_folders.all,
        regex: /sample-42/)

      # Should create the destination folder
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)

      # Renaming the entire directory tree doesn't require deleting the source ones
      expect(WebMock).not_to have_requested(:delete, /#{METIS_HOST}\/#{PROJECT}\/folder\/remove\/#{RELEASE_BUCKET}/)

      # Rename directory
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/rename\/#{RELEASE_BUCKET}/)
    end

    it 'when parent folder exists' do
      stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RESTRICT_BUCKET}\//)
        .to_return({
          status: 200,  # so the recursive method is called
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [],
            files: []
          }.to_json
        }).times(4).then
        .to_return({
          status: 422  # so the directory is created
        })
      stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\//)
        .to_return({
          status: 200,
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [{
              "folder_name": "premessa",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "folder_path": "assay/processed/sample-10/premessa"
            }],
            files: [{
              "file_name": "data.txt",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "file_path": "assay/processed/sample-10/data.txt"
            }]
          }.to_json}).then
        .to_return({
          status: 200,
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [],
            files: [{
              "file_name": "more-data.txt",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "file_path": "assay/processed/sample-10/premessa/more-data.txt"
            }]
          }.to_json})
      stub_create_folder({bucket: RESTRICT_BUCKET})
      stub_rename_file({bucket: RELEASE_BUCKET})
      stub_delete_folder({bucket: RELEASE_BUCKET})

      source_folders = Etna::Clients::Metis::Folders.new([{
        "folder_name": "sample-10",
        "bucket_name": RELEASE_BUCKET,
        "project_name": PROJECT,
        "folder_path": "assay/processed/sample-10"
      }])

      test_class.rename_folders_by_regex(
        source_bucket: RELEASE_BUCKET,
        project_name: PROJECT,
        dest_bucket: RESTRICT_BUCKET,
        source_folders: source_folders.all,
        regex: /sample-10/)

      # Should create the dest folder
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/folder\/create\/#{RESTRICT_BUCKET}/)

      # Remove the two source folders, sample-10 and premessa
      expect(WebMock).to have_requested(:delete, /#{METIS_HOST}\/#{PROJECT}\/folder\/remove\/#{RELEASE_BUCKET}/).times(2)

      # Rename both nested files
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/file\/rename\/#{RELEASE_BUCKET}/).times(2)
    end

    it 'when parent folder exists and some file exists in destination' do
      stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RESTRICT_BUCKET}\//)
        .to_return({
          status: 200,  # so the recursive method is called
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [],
            files: []
          }.to_json
        }).times(4).then
        .to_return({
          status: 200, # List an existing file
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [],
            files: [{
              "file_name": "data.txt",
              "bucket_name": RESTRICT_BUCKET,
              "project_name": PROJECT,
              "file_path": "assay/processed/sample-10/data.txt",
              "updated_at": "2000-02-01 00:00:00"
            }]
          }.to_json
        }).times(2).then
        .to_return({
          status: 422
        })
      stub_request(:get, /#{METIS_HOST}\/#{PROJECT}\/list\/#{RELEASE_BUCKET}\//)
        .to_return({
          status: 200,
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [{
              "folder_name": "premessa",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "folder_path": "assay/processed/sample-10/premessa"
            }],
            files: [{
              "file_name": "data.txt",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "file_path": "assay/processed/sample-10/data.txt",
              "updated_at": "2000-01-01 00:00:00"
            }]
          }.to_json}).then
        .to_return({
          status: 200,
          headers: {
            'Content-Type': 'application/json'
          },
          body: {
            folders: [],
            files: [{
              "file_name": "more-data.txt",
              "bucket_name": RELEASE_BUCKET,
              "project_name": PROJECT,
              "file_path": "assay/processed/sample-10/premessa/more-data.txt",
              "updated_at": "2000-01-01 00:00:00"
            }]
          }.to_json})
      stub_create_folder({bucket: RESTRICT_BUCKET})
      stub_rename_file({bucket: RELEASE_BUCKET})
      stub_delete_folder({bucket: RELEASE_BUCKET})
      stub_delete_file({bucket: RELEASE_BUCKET})

      source_folders = Etna::Clients::Metis::Folders.new([{
        "folder_name": "sample-10",
        "bucket_name": RELEASE_BUCKET,
        "project_name": PROJECT,
        "folder_path": "assay/processed/sample-10"
      }])

      test_class.rename_folders_by_regex(
        source_bucket: RELEASE_BUCKET,
        project_name: PROJECT,
        dest_bucket: RESTRICT_BUCKET,
        source_folders: source_folders.all,
        regex: /sample-10/)

      # Remove the two source folders, sample-10 and premessa
      expect(WebMock).to have_requested(:delete, /#{METIS_HOST}\/#{PROJECT}\/folder\/remove\/#{RELEASE_BUCKET}/).times(2)

      # Rename one, non-existing file
      expect(WebMock).to have_requested(:post, /#{METIS_HOST}\/#{PROJECT}\/file\/rename\/#{RELEASE_BUCKET}/)

      # Delete the existing file from source
      expect(WebMock).to have_requested(:delete, /#{METIS_HOST}\/#{PROJECT}\/file\/remove\/#{RELEASE_BUCKET}/)
    end

  end
end