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

    def rename_file(path, new_path)
      json_post("/athena/file/rename/files/#{path}", new_file_path: new_path)
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
      expect(json_body[:error]).to eq('Invalid path')
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

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File not found')

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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot rename over existing file')

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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot rename over existing folder')

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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('File is read-only')
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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Folder is read-only')
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file.folder).to be_nil
      expect(@wisdom_file).to be_has_data
    end

    it 'will not move a file to a non-existent folder' do
      token_header(:editor)
      rename_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid folder')
      @wisdom_file.refresh
      expect(@wisdom_file.file_path).to eq('wisdom.txt')
      expect(@wisdom_file.folder).to be_nil
      expect(@wisdom_file).to be_has_data
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

      # there is a new file
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
      expect(json_body[:error]).to eq("Invalid path for dest metis://athena/files/learn\nwisdom.txt")

      # the original is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
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

      # there is no new file
      expect(Metis::File.count).to eq(1)
    end

    it 'refuses to copy a non-existent file' do
      # we attempt to rename a file that does not exist
      token_header(:editor)
      copy_file('folly.txt', 'learn-folly.txt')

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File metis://athena/files/folly.txt not found')

      # the actual file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
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

      # the file we tried to rename is untouched
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')

      # we can still see the data
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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot copy over existing folder metis://athena/files/learn-wisdom.txt')

      # the file we tried to copy is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we tried to copy is untouched
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
      copy_file('wisdom.txt', 'contents/wisdom.txt')

      expect(last_response.status).to eq(200)

      # the original is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # a new file exists
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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('contents folder is read-only')

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
      expect(json_body[:error]).to eq('Invalid folder')

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

      # there is a new file
      expect(Metis::File.count).to eq(2)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(sundry_bucket)
    end

    it 'refuses to copy without bucket permissions' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'administrator', owner: 'metis' )
      copy_file('wisdom.txt', 'learn-wisdom.txt', new_bucket_name: 'sundry')
      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot access bucket')

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
      expect(json_body[:error]).to eq('Cannot access bucket')

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
      json_post("/athena/file/bulk_copy", {
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

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: "metis://athena/files/learn\nwisdom.txt"
      }])

      @wisdom_file.refresh
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('All revisions require valid "dest" parameter')

      # the original is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      orig_helmet_file = Metis::File.last
      expect(orig_helmet_file.file_name).to eq('helmet.jpg')
    end

    it 'refuses to copy a file without permissions' do
      # the user is a viewer, not an editor
      token_header(:viewer)

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
      # we attempt to rename a file that does not exist
      token_header(:editor)

      bulk_copy([{
        source: 'metis://athena/files/folly.txt',
        dest: 'metis://athena/files/learn-folly.txt'
      }])

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File metis://athena/files/folly.txt not found')

      # the actual file is untouched
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

      # the file we tried to rename is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we tried to rename now points to the new location
      new_link_file = Metis::File.last
      expect(new_link_file.file_name).to eq('learn-wisdom.txt')
      expect(new_link_file.data_block.location).to eq(@helmet_file.data_block.location)

      # we can still see the data
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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot copy over existing folder metis://athena/files/learn-wisdom.txt')

      # the file we tried to copy is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')

      # the file we tried to copy is untouched
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

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('contents folder is read-only')

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

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/contents/learn-wisdom.txt'
      }])

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid folder')

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

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }])

      stubs.add_file('athena', 'sundry', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is a new file
      expect(Metis::File.count).to eq(3)
      new_wisdom_file = Metis::File.last
      expect(new_wisdom_file.file_name).to eq('learn-wisdom.txt')
      expect(new_wisdom_file).to be_has_data
      expect(new_wisdom_file.data_block).to eq(@wisdom_file.data_block)

      expect(new_wisdom_file.bucket).to eq(sundry_bucket)
    end

    it 'refuses to copy without bucket permissions' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'administrator', owner: 'metis' )

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }])

      stubs.add_file('athena', 'files', 'learn-wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('User does not have access to all dest buckets')

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
      expect(Metis::File.count).to eq(2)
    end

    it 'copies to a bucket with a different application owner with matching signature' do
      token_header(:editor)
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'vulcan' )

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }], hmac_params(id: 'vulcan', signature: 'valid'))

      stubs.add_file('athena', 'sundry', 'learn-wisdom.txt')

      expect(last_response.status).to eq(200)

      # the old file is untouched
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data

      # there is no new file
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

      bulk_copy([{
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/sundry/learn-wisdom.txt'
      }])

      stubs.add_file('athena', 'sundry', 'learn-wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('User does not have access to all dest buckets')

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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('contents folder is read-only')

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

    it 'creates multiple links in a single update' do
      token_header(:editor)

      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

      # there is a new file
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
      expect(Metis::File.count).to eq(5)

      # the data is not destroyed
      expect(::File.exists?(location)).to be_truthy
      orig_wisdom_file = Metis::File.first
      expect(orig_wisdom_file.file_name).to eq('wisdom.txt')
      third_link_file = Metis::File.last
      expect(third_link_file.file_name).to eq('build-helmet.jpg')
      expect(third_link_file.bucket.name).to eq('sundry')
    end
  end
end
