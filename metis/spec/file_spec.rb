require 'digest'

describe FileController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    default_bucket('athena')

    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    stubs.clear
  end

  context '#remove' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')

      @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)
      stubs.create_folder('athena', 'files', 'blueprints/helmet')

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @helmet_folder)
      stubs.create_file('athena', 'files', 'blueprints/helmet/helmet.jpg', HELMET)
    end

    def remove_file(path)
      delete("athena/file/remove/files/#{path}")
    end

    it 'removes a file' do
      token_header(:editor)
      location = @helmet_file.data_block.location
      remove_file('blueprints/helmet/helmet.jpg')

      expect(last_response.status).to eq(200)
      expect(Metis::File.count).to eq(1)

      # the data is not destroyed
      expect(::File.exists?(location)).to be_truthy

      location = @wisdom_file.data_block.location
      remove_file('wisdom.txt')

      # the data is not destroyed
      expect(::File.exists?(location)).to be_truthy
      expect(last_response.status).to eq(200)
      expect(Metis::File.count).to eq(0)
    end

    it 'refuses to remove a file without permissions' do
      # we attempt to remove a file though we are a mere viewer
      token_header(:viewer)
      remove_file('wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(@wisdom_file).to be_has_data
      expect(Metis::File.count).to eq(2)
    end

    it 'refuses to remove a non-existent file' do
      # we attempt to remove a file that does not exist
      token_header(:editor)
      remove_file('folly.txt')

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File not found')
    end

    it 'refuses to remove a read-only file' do
      @wisdom_file.read_only = true
      @wisdom_file.save
      @wisdom_file.refresh

      token_header(:editor)
      remove_file('wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('File is read-only')
      expect(@wisdom_file).to be_has_data
      expect(Metis::File.all).to include(@wisdom_file)
    end

    it 'refuses to remove a read-only file even for an admin' do
      @wisdom_file.read_only = true
      @wisdom_file.save
      @wisdom_file.refresh

      token_header(:admin)
      remove_file('wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('File is read-only')
      expect(@wisdom_file).to be_has_data
      expect(Metis::File.all).to include(@wisdom_file)
    end
  end

  context '#protect' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
    end

    def protect_file path, params={}
      json_post("/athena/file/protect/files/#{path}", params)
    end

    it 'protects a file' do
      token_header(:admin)
      protect_file('wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(200)
      expect(json_body[:files].first[:read_only]).to eq(true)
      expect(@wisdom_file).to be_read_only
    end

    it 'refuses to protect a file without permissions' do
      token_header(:editor)
      protect_file('wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(403)
      expect(@wisdom_file).not_to be_read_only
    end

    it 'refuses to protect a non-existent file' do
      # we attempt to protect a file that does not exist
      token_header(:admin)
      protect_file('folly.txt')

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File not found')

      # the actual file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file).not_to be_read_only
    end

    it 'refuses to protect a read-only file' do
      @wisdom_file.read_only = true
      @wisdom_file.save
      @wisdom_file.refresh

      token_header(:admin)
      protect_file('wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('File is already read-only')
      @wisdom_file.refresh
      expect(@wisdom_file).to be_read_only
    end
  end

  context '#unprotect' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, read_only: true)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
      expect(@wisdom_file).to be_read_only
    end

    def unprotect_file path, params={}
      json_post("/athena/file/unprotect/files/#{path}",params)
    end

    it 'unprotects a file' do
      token_header(:admin)
      unprotect_file('wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(200)
      expect(@wisdom_file).not_to be_read_only
    end

    it 'refuses to unprotect a file without permissions' do
      token_header(:editor)
      unprotect_file('wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(403)
      expect(@wisdom_file).to be_read_only
    end

    it 'refuses to unprotect a non-existent file' do
      # we attempt to unprotect a file that does not exist
      token_header(:admin)
      unprotect_file('folly.txt')

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File not found')

      # the actual file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file).to be_read_only
    end

    it 'refuses to unprotect a writeable file' do
      @wisdom_file.read_only = false
      @wisdom_file.save

      token_header(:admin)
      unprotect_file('wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('File is not protected')
      @wisdom_file.refresh
      expect(@wisdom_file).not_to be_read_only
    end
  end

  context '#rename' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
    end

    def rename_file(path, new_path, params={})
      json_post("/athena/file/rename/files/#{path}", {
        new_file_path: new_path
      }.merge(params))
    end

    it 'renames a file' do
      token_header(:editor)
      rename_file('wisdom.txt', 'learn-wisdom.txt')
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(200)
      expect(@wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'refuses to rename a file to an invalid name' do
      token_header(:editor)
      rename_file('wisdom.txt', "learn\nwisdom.txt")

      @wisdom_file.refresh
      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq(["Invalid path: \"metis://athena/files/learn\nwisdom.txt\""])
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'refuses to rename a file without permissions' do
      # the user is a viewer, not an editor
      token_header(:viewer)
      rename_file('wisdom.txt', 'learn-wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(403)
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'refuses to rename a non-existent file' do
      # we attempt to rename a file that does not exist
      token_header(:editor)
      rename_file('folly.txt', 'learn-folly.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq(["File \"metis://athena/files/folly.txt\" not found"])

      # the actual file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'refuses to rename over an existing file' do
      learn_wisdom_file = create_file('athena', 'learn-wisdom.txt', WISDOM*2)
      stubs.create_file('athena', 'files', 'learn-wisdom.txt', WISDOM*2)

      token_header(:editor)
      rename_file('wisdom.txt', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq(
        ["File \"metis://athena/files/learn-wisdom.txt\" already exists"])

      # the file we tried to rename is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we tried to rename is untouched
      learn_wisdom_file.refresh
      expect(Metis::File.last).to eq(learn_wisdom_file)
      expect(learn_wisdom_file.file_name).to eq('learn-wisdom.txt')

      # we can still see the data
      expect(@wisdom_file).to be_has_data
      expect(learn_wisdom_file).to be_has_data
      expect(::File.read(@wisdom_file.data_block.location)).to eq(WISDOM)
      expect(::File.read(learn_wisdom_file.data_block.location)).to eq(WISDOM*2)
    end

    it 'refuses to rename over an existing folder' do
      learn_wisdom_folder = create_folder('athena', 'learn-wisdom.txt')
      stubs.create_folder('athena', 'files', 'learn-wisdom.txt')

      token_header(:editor)
      rename_file('wisdom.txt', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq(
        ["Cannot copy over existing folder: \"metis://athena/files/learn-wisdom.txt\""])

      # the file we tried to rename is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we tried to rename is untouched
      learn_wisdom_folder.refresh
      expect(Metis::Folder.last).to eq(learn_wisdom_folder)
      expect(learn_wisdom_folder.folder_name).to eq('learn-wisdom.txt')

      # we can still see the data
      expect(@wisdom_file).to be_has_data
      expect(::File.read(@wisdom_file.data_block.location)).to eq(WISDOM)
    end

    it 'refuses to rename a read-only file' do
      @wisdom_file.read_only = true
      @wisdom_file.save

      token_header(:editor)
      rename_file('wisdom.txt', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq(
        ["File \"metis://athena/files/wisdom.txt\" is read-only"])
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'can move a file to a new folder' do
      contents_folder = create_folder('athena', 'contents')
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)
      rename_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(200)
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('contents/wisdom.txt')
      expect(@wisdom_file.folder).to eq(contents_folder)
      expect(@wisdom_file).to be_has_data
    end

    it 'will not move a file to a read-only folder' do
      contents_folder = create_folder('athena', 'contents', read_only: true)
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)
      rename_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq(["Folder \"contents\" is read-only"])
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file.folder).to be_nil
      expect(@wisdom_file).to be_has_data
    end

    it 'will not move a file to a non-existent folder' do
      token_header(:editor)
      rename_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid folder: "contents"')
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file.folder).to be_nil
      expect(@wisdom_file).to be_has_data
    end

    it 'renames a file into a new bucket and folder' do
      token_header(:editor)
      sundry_bucket = create(:bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )
      sundry_wisdom_folder = create_folder('athena', 'sundry-wisdom', bucket: sundry_bucket)
      stubs.create_folder('athena', 'sundry', 'sundry-wisdom')

      rename_file('wisdom.txt', 'sundry-wisdom/updated-wisdom.txt', {
        new_bucket_name: 'sundry'
      })

      expect(last_response.status).to eq(200)

      # the old file is updated
      expect(Metis::File.count).to eq(1)
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('updated-wisdom.txt')
      expect(@wisdom_file).to be_has_data
      expect(@wisdom_file.bucket).to eq(sundry_bucket)
      expect(@wisdom_file.folder).to eq(sundry_wisdom_folder)
    end

    it 'renames a file into a new bucket and no folder' do
      token_header(:editor)
      sundry_bucket = create(:bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

      rename_file('wisdom.txt', 'updated-wisdom.txt', {
        new_bucket_name: 'sundry'
      })

      expect(last_response.status).to eq(200)

      # the old file is updated
      expect(Metis::File.count).to eq(1)
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('updated-wisdom.txt')
      expect(@wisdom_file).to be_has_data
      expect(@wisdom_file.bucket).to eq(sundry_bucket)
      expect(@wisdom_file.folder).to eq(nil)
    end
  end

  context '#copy' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
    end

    def copy_file(path, new_path, params={})
      json_post("/athena/file/copy/files/#{path}", {
        new_file_path: new_path
      }.merge(params))
    end

    it 'copies a file' do
      token_header(:editor)
      copy_file('wisdom.txt', 'learn-wisdom.txt')
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file that points to the same data block
      expect(Metis::File.count).to eq(2)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)
    end

    it 'refuses to copy a file to an invalid name' do
      token_header(:editor)
      copy_file('wisdom.txt', "learn\nwisdom.txt")

      @wisdom_file.refresh
      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        "Invalid path: \"metis://athena/files/learn\nwisdom.txt\""
      )

      # the original is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file created
      expect(Metis::File.count).to eq(1)
    end

    it 'refuses to copy a file without permissions' do
      # the user is a viewer, not an editor
      token_header(:viewer)
      copy_file('wisdom.txt', 'learn-wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(403)

      # the original is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file created
      expect(Metis::File.count).to eq(1)
    end

    it 'refuses to copy a non-existent file' do
      # we attempt to rename a file that does not exist
      token_header(:editor)
      copy_file('folly.txt', 'learn-folly.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        "File \"metis://athena/files/folly.txt\" not found"
      )

      # the actual file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file created
      expect(Metis::File.count).to eq(1)
    end

    it 'can replace an existing copy' do
      learn_wisdom_file = create_file('athena', 'learn-wisdom.txt', WISDOM*2)
      stubs.create_file('athena', 'files', 'learn-wisdom.txt', WISDOM*2)

      token_header(:editor)
      copy_file('wisdom.txt', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the file we tried to rename is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we updated keeps the old filename
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')

      # The file we updated points to the new data block
      expect(@wisdom_file).to be_has_data
      expect(new_wisdom_file).to be_has_data
      expect(::File.read(@wisdom_file.data_block.location)).to eq(WISDOM)
      expect(new_wisdom_file.data_block.location).to eq(@wisdom_file.data_block.location)
    end

    it 'refuses to copy over an existing folder' do
      learn_wisdom_folder = create_folder('athena', 'learn-wisdom.txt')
      stubs.create_folder('athena', 'files', 'learn-wisdom.txt')

      token_header(:editor)
      copy_file('wisdom.txt', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        "Cannot copy over existing folder: \"metis://athena/files/learn-wisdom.txt\""
      )

      # the file we tried to copy is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the folder we tried to overwrite is untouched
      learn_wisdom_folder.refresh
      expect(Metis::Folder.last).to eq(learn_wisdom_folder)
      expect(learn_wisdom_folder.folder_name).to eq('learn-wisdom.txt')

      # we can still see the original file's data
      expect(@wisdom_file).to be_has_data
      expect(::File.read(@wisdom_file.data_block.location)).to eq(WISDOM)
    end

    it 'can copy a file to a new folder' do
      contents_folder = create_folder('athena', 'contents')
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)
      copy_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(200)

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # a new file exists in a new folder
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_path).to eq('contents/wisdom.txt')
      expect(new_wisdom_file.folder).to eq(contents_folder)
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)
    end

    it 'will not copy a file to a read-only folder' do
      contents_folder = create_folder('athena', 'contents', read_only: true)
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)
      copy_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        "Folder \"contents\" is read-only"
      )

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(1)
    end

    it 'will not copy a file to a non-existent folder' do
      token_header(:editor)
      copy_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("Invalid folder: \"contents\"")

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(1)
    end

    it 'copies a file to a new bucket' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )
      copy_file('wisdom.txt', 'learn-wisdom.txt', new_bucket_name: 'sundry')
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file in the new bucket
      expect(Metis::File.count).to eq(2)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(sundry_bucket)
    end

    it 'copies a file to a different project' do
      backup_bucket = default_bucket('backup')
      token_header(:editor)

      expect(Metis::File.count).to eq(1)

      copy_file('wisdom.txt', 'metis://backup/files/learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file in the new bucket
      expect(Metis::File.count).to eq(2)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(backup_bucket)
      expect(new_wisdom_file.project_name).to eq('backup')
    end

    it 'throws exception if no permissions on dest project' do
      backup_bucket = default_bucket('backup')

      token_header(:editor)

      expect(Metis::File.count).to eq(1)

      copy_file('wisdom.txt', 'metis://second_backup/files/learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq(
        "Invalid bucket on project second_backup: \"files\""
      )

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(1)
    end

    it 'refuses to copy without bucket permissions' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'administrator', owner: 'metis' )
      copy_file('wisdom.txt', 'learn-wisdom.txt', new_bucket_name: 'sundry')
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq(
        "Cannot access bucket: \"sundry\""
      )

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(1)
    end

    it 'copies to a bucket with a different application owner with matching signature' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'vulcan' )
      copy_file('wisdom.txt', 'learn-wisdom.txt', {new_bucket_name: 'sundry'}.merge(
        hmac_params(id: 'vulcan', signature: 'valid')))
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
    end

    it 'refuses to copy to a bucket with a different owner without a signature' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'vulcan' )
      copy_file('wisdom.txt', 'learn-wisdom.txt',
        hmac_params.merge(signature: 'valid', new_bucket_name: 'sundry')
      )
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq(
        "Cannot access bucket: \"sundry\""
      )

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(1)
    end
  end

  context '#bulk_copy' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')

      @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)
      stubs.create_folder('athena', 'files', 'blueprints/helmet')

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @helmet_folder)
      stubs.create_file('athena', 'files', 'blueprints/helmet/helmet.jpg', HELMET)
    end

    def bulk_copy(revisions=[], params={})
      json_post("/athena/files/copy", {
        revisions: revisions
      }.merge(params))
    end

    def copy_file(path, new_path, params={})
      json_post("/athena/file/copy/files/#{path}", {
        new_file_path: new_path
      }.merge(params))
    end

    it 'copies a file' do
      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(200)

      # the old file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file
      expect(Metis::File.count).to eq(3)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)
    end

    it 'refuses to copy a file to an invalid name' do
      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: "metis://athena/files/learn\nwisdom.txt"
      }])

      @wisdom_file.refresh
      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/files/learn\nwisdom.txt",
         "errors": ["Invalid path: \"metis://athena/files/learn\nwisdom.txt\""],
         "source": "metis://athena/files/wisdom.txt"}
      )

      # the original is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'refuses to copy a file without permissions' do
      # the user is a viewer, not an editor
      token_header(:viewer)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/learn-wisdom.txt'
      }])

      @wisdom_file.refresh
      expect(last_response.status).to eq(403)

      # the original is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'refuses to copy a non-existent file' do
      # we attempt to copy a file that does not exist
      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/folly.txt',
        dest: 'metis://athena/files/learn-folly.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/files/learn-folly.txt",
         "errors": ["File \"metis://athena/files/folly.txt\" not found"],
         "source": "metis://athena/files/folly.txt"}
      )

      # no new file is created
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'can replace an existing copy' do
      learn_wisdom_file = create_file('athena', 'learn-wisdom.txt', WISDOM*2)
      stubs.create_file('athena', 'files', 'learn-wisdom.txt', WISDOM*2)

      expect(Metis::File.count).to eq(3)

      token_header(:editor)
      bulk_copy([{
        source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
        dest: 'metis://athena/files/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(200)

      expect(Metis::File.count).to eq(3)

      # the file we copied is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we updated now points to the new location
      new_link_file = Metis::File.last
      expect(new_link_file.file_name).to eq('learn-wisdom.txt')
      expect(new_link_file.data_block.location).to eq(@helmet_file.data_block.location)

      # we can still see the original data
      expect(@wisdom_file).to be_has_data
      expect(new_link_file).to be_has_data
      expect(::File.read(@wisdom_file.data_block.location)).to eq(WISDOM)
    end

    it 'refuses to copy over an existing folder' do
      learn_wisdom_folder = create_folder('athena', 'learn-wisdom.txt')
      stubs.create_folder('athena', 'files', 'learn-wisdom.txt')

      token_header(:editor)
      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/files/learn-wisdom.txt",
        "errors": ["Cannot copy over existing folder: \"metis://athena/files/learn-wisdom.txt\""],
        "source": "metis://athena/files/wisdom.txt"}
      )

      # the file we tried to copy is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the folder we tried to overwrite is untouched
      learn_wisdom_folder.refresh
      expect(Metis::Folder.last).to eq(learn_wisdom_folder)
      expect(learn_wisdom_folder.folder_name).to eq('learn-wisdom.txt')

      # we can still see the data
      expect(@wisdom_file).to be_has_data
      expect(::File.read(@wisdom_file.data_block.location)).to eq(WISDOM)
    end

    it 'can copy a file to a new folder' do
      contents_folder = create_folder('athena', 'contents')
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(200)

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # a new file exists
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_path).to eq('contents/learn-wisdom.txt')
      expect(new_wisdom_file.folder).to eq(contents_folder)
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)
    end

    it 'will not copy a file to a read-only folder' do
      contents_folder = create_folder('athena', 'contents', read_only: true)
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/files/contents/learn-wisdom.txt",
         "errors": ["Folder \"contents\" is read-only"],
         "source": "metis://athena/files/wisdom.txt"}
      )

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'will not copy a file to a non-existent folder' do
      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/files/contents/learn-wisdom.txt",
         "errors": ["Invalid folder: \"contents\""],
         "source": "metis://athena/files/wisdom.txt"}
      )

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'copies a file to a new bucket' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }])

      stubs.add_file('athena', 'sundry', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file in the new bucket
      expect(Metis::File.count).to eq(3)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(sundry_bucket)
    end

    it 'copies a file to a new project' do
      backup_bucket = default_bucket('backup')

      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://backup/files/learn-wisdom.txt'
      }])

      stubs.add_file('backup', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file in the new bucket
      expect(Metis::File.count).to eq(3)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(backup_bucket)
      expect(new_wisdom_file.project_name).to eq('backup')
    end

    it 'throws exception if no permissions on dest project' do
      backup_bucket = default_bucket('backup')

      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://second_backup/files/learn-wisdom.txt'
      }])

      stubs.add_file('backup', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors]).to eq([{
        dest: "metis://second_backup/files/learn-wisdom.txt",
        source: "metis://athena/files/wisdom.txt",
        errors: ["Invalid bucket \"files\" in project \"second_backup\". Check the bucket name and your permissions."]}])

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file in the new bucket
      expect(Metis::File.count).to eq(2)
    end

    it 'refuses to copy without bucket permissions' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'administrator', owner: 'metis' )

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }])

      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/sundry/learn-wisdom.txt",
         "errors": ["Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."],
         "source": "metis://athena/files/wisdom.txt"}
      )

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
    end

    it 'copies to a bucket with a different application owner with matching signature' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'vulcan' )

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }], hmac_params(id: 'vulcan', signature: 'valid'))

      stubs.add_file('athena', 'sundry', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file in the new bucket
      expect(Metis::File.count).to eq(3)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(sundry_bucket)
    end

    it 'refuses to copy to a bucket with a different owner without a signature' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'vulcan' )

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }])

      stubs.add_file('athena', 'sundry', 'learn-wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/sundry/learn-wisdom.txt",
         "errors": ["Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."],
         "source": "metis://athena/files/wisdom.txt"}
      )

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
    end

    it 'refuses to update if a single revision out of multiple is invalid' do
      contents_folder = create_folder('athena', 'contents', read_only: true)
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/learn-wisdom.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/learn-wisdom2.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(1)
      expect(json_body[:errors][0]).to eq(
        {"dest": "metis://athena/files/contents/learn-wisdom.txt",
         "errors": ["Folder \"contents\" is read-only"],
         "source": "metis://athena/files/wisdom.txt"}
      )

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'creates multiple copies in a single update' do
      token_header(:editor)

      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

      expect(Metis::File.count).to eq(2)

      location = @wisdom_file.data_block.location
      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/learn-wisdom2.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom3.txt'
      }, {
        source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
        dest: 'metis://athena/sundry/build-helmet.jpg'
      }])

      expect(last_response.status).to eq(200)

      # the data is not destroyed
      expect(::File.exists?(location)).to be_truthy

      # There are some new files
      expect(Metis::File.count).to eq(5)
      third_link_file = Metis::File.find(file_name: 'build-helmet.jpg')
      expect(third_link_file.file_name).to eq('build-helmet.jpg')
      expect(third_link_file.bucket.name).to eq('sundry')
    end

    it 'can swap data blocks during copy (i.e. caches pre-copy data blocks)' do
      # @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      # stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      # @blueprints_folder = create_folder('athena', 'blueprints')
      # stubs.create_folder('athena', 'files', 'blueprints')

      # @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)
      # stubs.create_folder('athena', 'files', 'blueprints/helmet')

      # @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @helmet_folder)
      # stubs.create_file('athena', 'files', 'blueprints/helmet/helmet.jpg', HELMET)
      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      original_wisdom_block = @wisdom_file.data_block
      original_helmet_block = @helmet_file.data_block
      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/blueprints/helmet/helmet.jpg'
      }, {
        source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
        dest: 'metis://athena/files/wisdom.txt'
      }])

      expect(last_response.status).to eq(200)

      # the data is not destroyed
      expect(::File.exists?(original_wisdom_block.location)).to be_truthy
      expect(::File.exists?(original_helmet_block.location)).to be_truthy

      # There are some new files
      expect(Metis::File.count).to eq(2)
      @wisdom_file.refresh
      @helmet_file.refresh

      expect(@wisdom_file.data_block).to eq(original_helmet_block)
      expect(@helmet_file.data_block).to eq(original_wisdom_block)
    end

    it 'receives all errors back for multiple invalid revisions' do
      contents_folder = create_folder('athena', 'contents', read_only: true)
      stubs.create_folder('athena', 'files', 'contents')

      token_header(:editor)

      expect(Metis::File.count).to eq(2)

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: nil
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom2.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:errors].length).to eq(3)
      expect(json_body[:errors]).to eq(
        [{"dest": "metis://athena/files/contents/learn-wisdom.txt",
          "errors": ["Folder \"contents\" is read-only"],
          "source": "metis://athena/files/wisdom.txt"},
         {"dest": nil,
          "errors": ["Invalid path: \"\""],
          "source": "metis://athena/files/wisdom.txt"},
         {"dest": "metis://athena/sundry/learn-wisdom2.txt",
          "errors": ["Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."],
          "source": "metis://athena/files/wisdom.txt"}]
      )

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is are no new files
      expect(Metis::File.count).to eq(2)
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'can copy files from multiple buckets to multiple buckets in single operation' do
      contents_folder = create_folder('athena', 'contents')
      stubs.create_folder('athena', 'files', 'contents')
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )
      magma_bucket = create( :bucket, project_name: 'athena', name: 'magma', access: 'administrator', owner: 'magma' )

      shiny_helmet_file = create_file('athena', 'shiny-helmet.jpg', SHINY_HELMET, bucket: sundry_bucket)
      stubs.create_file('athena', 'sundry', 'shiny-helmet.jpg', SHINY_HELMET)

      token_header(:editor)

      expect(Metis::File.count).to eq(3)

      location = @wisdom_file.data_block.location
      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/magma/learn-wisdom2.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom3.txt'
      }, {
        source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
        dest: 'metis://athena/sundry/build-helmet.jpg'
      }, {
        source: 'metis://athena/sundry/shiny-helmet.jpg',
        dest: 'metis://athena/magma/polish-helmet.jpg'
      }], hmac_params(id: 'magma', signature: 'valid'))

      expect(last_response.status).to eq(200)

      # the data is not destroyed
      expect(::File.exists?(location)).to be_truthy

      # There are new files
      expect(Metis::File.count).to eq(8)
      orig_wisdom_file = Metis::File.find(file_name: 'wisdom.txt')
      expect(orig_wisdom_file.data_block.location).to eq(location)
      learn_wisdom_2_file = Metis::File.find(file_name: 'learn-wisdom2.txt')
      expect(learn_wisdom_2_file.bucket.name).to eq('magma')
      build_helmet_file = Metis::File.find(file_name: 'build-helmet.jpg')
      expect(build_helmet_file.bucket.name).to eq('sundry')
    end

    it 'refuses to copy from multiple buckets to multiple buckets if does not have permissions for one' do
      contents_folder = create_folder('athena', 'contents')
      stubs.create_folder('athena', 'files', 'contents')
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'administrator', owner: 'metis' )
      magma_bucket = create( :bucket, project_name: 'athena', name: 'magma', access: 'administrator', owner: 'magma' )

      shiny_helmet_file = create_file('athena', 'shiny-helmet.jpg', SHINY_HELMET, bucket: sundry_bucket)
      stubs.create_file('athena', 'sundry', 'shiny-helmet.jpg', SHINY_HELMET)

      token_header(:editor)

      expect(Metis::File.count).to eq(3)

      location = @wisdom_file.data_block.location
      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/magma/learn-wisdom2.txt'
      }, {
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom3.txt'
      }, {
        source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
        dest: 'metis://athena/sundry/build-helmet.jpg'
      }, {
        source: 'metis://athena/sundry/shiny-helmet.jpg',
        dest: 'metis://athena/magma/polish-helmet.jpg'
      }], hmac_params(id: 'magma', signature: 'valid'))

      expect(last_response.status).to eq(422)

      expect(json_body[:errors].length).to eq(3)
      expect(json_body[:errors]).to eq(
        [{"dest": "metis://athena/sundry/learn-wisdom3.txt",
         "errors": ["Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."],
         "source": "metis://athena/files/wisdom.txt"},
        {"dest": "metis://athena/sundry/build-helmet.jpg",
         "errors": ["Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."],
         "source": "metis://athena/files/blueprints/helmet/helmet.jpg"},
        {"dest": "metis://athena/magma/polish-helmet.jpg",
         "errors": ["Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."],
         "source": "metis://athena/sundry/shiny-helmet.jpg"}])

      # there are no new files
      expect(Metis::File.count).to eq(3)

      # the data is not destroyed
      expect(::File.exists?(location)).to be_truthy
    end
  end

  context '#update_bucket!' do
    before(:each) do
      @backup_files_bucket = create( :bucket, project_name: 'athena', name: 'backup_files', owner: 'metis', access: 'viewer')
      stubs.create_bucket('athena', 'backup_files')

      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')

      @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)
      stubs.create_folder('athena', 'files', 'blueprints/helmet')

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @helmet_folder)
      stubs.create_file('athena', 'files', 'blueprints/helmet/helmet.jpg', HELMET)
    end

    it 'moves a file to the new bucket when does not have a parent folder' do
      expect(@wisdom_file.bucket).not_to eq(@backup_files_bucket)

      @wisdom_file.update_bucket!(@backup_files_bucket)

      expect(@wisdom_file.bucket).to eq(@backup_files_bucket)
    end

    it 'moves a file to the new bucket when does have a parent folder' do
      backup_blueprints_folder = create_folder('athena', 'blueprints', bucket: @backup_files_bucket)
      stubs.create_folder('athena', 'backup_files', 'blueprints')

      expect(@helmet_file.bucket).not_to eq(@backup_files_bucket)
      @helmet_file.folder = backup_blueprints_folder

      @helmet_file.update_bucket!(@backup_files_bucket)

      expect(@helmet_file.bucket).to eq(@backup_files_bucket)
    end

    it 'cannot move a file to a different bucket than its folder\'s bucket' do
      expect(@helmet_file.bucket).not_to eq(@backup_files_bucket)

      expect {
        @helmet_file.update_bucket!(@backup_files_bucket)
      }.to raise_error(StandardError)

      expect(@helmet_file.bucket).not_to eq(@backup_files_bucket)
    end
  end

  context '#update_bucket_and_rename!' do
    before(:each) do
      @backup_files_bucket = create( :bucket, project_name: 'athena', name: 'backup_files', owner: 'metis', access: 'viewer')
      stubs.create_bucket('athena', 'backup_files')

      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')

      @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)
      stubs.create_folder('athena', 'files', 'blueprints/helmet')

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @helmet_folder)
      stubs.create_file('athena', 'files', 'blueprints/helmet/helmet.jpg', HELMET)
    end

    it 'moves a file to the new bucket when does not have a parent folder' do
      expect(@wisdom_file.bucket).not_to eq(@backup_files_bucket)

      @wisdom_file.update_bucket_and_rename!(nil, 'additional_wisdom.txt', @backup_files_bucket)

      expect(@wisdom_file.bucket).to eq(@backup_files_bucket)
      expect(@wisdom_file.file_name).to eq('additional_wisdom.txt')
    end

    it 'moves a file to the new bucket when does have a parent folder' do
      backup_blueprints_folder = create_folder('athena', 'blueprints', bucket: @backup_files_bucket)
      stubs.create_folder('athena', 'backup_files', 'blueprints')

      expect(@helmet_file.bucket).not_to eq(@backup_files_bucket)

      @helmet_file.update_bucket_and_rename!(backup_blueprints_folder, 'backup_helmet.jpg', @backup_files_bucket)

      expect(@helmet_file.bucket).to eq(@backup_files_bucket)
      expect(@helmet_file.folder).to eq(backup_blueprints_folder)
      expect(@helmet_file.file_name).to eq('backup_helmet.jpg')
    end

    it 'cannot move a file to a different bucket than its folder\'s bucket' do
      expect(@helmet_file.bucket).not_to eq(@backup_files_bucket)

      expect {
        @helmet_file.update_bucket_and_rename!(@blueprints_folder, 'backup_helmet.jpg', @backup_files_bucket)
      }.to raise_error(StandardError)

      expect(@helmet_file.bucket).not_to eq(@backup_files_bucket)
      expect(@helmet_file.folder).not_to eq(@blueprints_folder)
      expect(@helmet_file.file_name).not_to eq('backup_helmet.jpg')
    end
  end

  context '#touch' do
    before(:each) do
      @creation_time = DateTime.now - 100
      Timecop.freeze(@creation_time)
      @blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @blueprints_folder)
      stubs.create_file('athena', 'files', 'blueprints/helmet.jpg', HELMET)
    end

    after(:each) do
      Timecop.return
    end

    def touch_file path
      get(
        "/athena/file/touch/files/#{path}"
      )
    end

    it 'updates a file timestamp' do
      expect(@helmet_file.updated_at.iso8601).to eq(@creation_time.iso8601)

      @update_time = DateTime.now
      Timecop.freeze(@update_time)

      token_header(:editor)
      touch_file('blueprints/helmet.jpg')

      @helmet_file.refresh
      expect(last_response.status).to eq(200)
      expect(@helmet_file.updated_at.iso8601).to eq(@update_time.iso8601)
    end

    it 'throws exception if file is read only' do
      expect(@helmet_file.updated_at.iso8601).to eq(@creation_time.iso8601)
      @helmet_file.read_only = true
      @helmet_file.save
      @helmet_file.refresh

      @update_time = DateTime.now
      Timecop.freeze(@update_time)

      token_header(:editor)
      touch_file('blueprints/helmet.jpg')

      @helmet_file.refresh
      expect(last_response.status).to eq(403)
      expect(@helmet_file.updated_at.iso8601).to eq(@creation_time.iso8601)
      
    end
  end
end
