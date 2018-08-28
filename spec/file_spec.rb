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
    clear_stubs
  end

  context '#remove' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stub_file('wisdom.txt', WISDOM, :athena)

      @blueprints_folder = create_folder('athena', 'blueprints')
      stub_folder('blueprints', 'athena')

      @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)
      stub_folder('blueprints/helmet', 'athena')

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET, folder: @helmet_folder)
      stub_file('blueprints/helmet/helmet.jpg', HELMET, :athena)
    end

    def remove_file(path)
      delete("athena/remove_file/files/#{path}")
    end


    it 'removes a file' do
      token_header(:editor)
      location = @helmet_file.location
      remove_file('blueprints/helmet/helmet.jpg')

      expect(last_response.status).to eq(200)
      expect(Metis::File.count).to eq(1)
      expect(::File.exists?(location)).to be_falsy

      location = @wisdom_file.location
      remove_file('wisdom.txt')

      expect(::File.exists?(location)).to be_falsy
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

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Cannot remove file')
      expect(@wisdom_file).to be_has_data
      expect(Metis::File.all).to include(@wisdom_file)
    end

    it 'refuses to remove a read-only file even for an admin' do
      @wisdom_file.read_only = true
      @wisdom_file.save
      @wisdom_file.refresh

      token_header(:admin)
      remove_file('wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Cannot remove file')
      expect(@wisdom_file).to be_has_data
      expect(Metis::File.all).to include(@wisdom_file)
    end
  end

  context '#protect' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stub_file('wisdom.txt', WISDOM, :athena)
    end

    def protect_file path, params={}
      json_post("athena/protect_file/files/#{path}", params)
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

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('File is read-only')
      @wisdom_file.refresh
      expect(@wisdom_file).to be_read_only
    end
  end

  context '#unprotect' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, read_only: true)
      stub_file('wisdom.txt', WISDOM, :athena)
      expect(@wisdom_file).to be_read_only
    end

    def unprotect_file path, params={}
      json_post("athena/unprotect_file/files/#{path}",params)
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
      stub_file('wisdom.txt', WISDOM, :athena)
    end

    def rename_file(path, new_path)
      json_post("athena/rename_file/files/#{path}", new_file_path: new_path)
    end

    it 'renames a file' do
      token_header(:editor)
      rename_file('wisdom.txt', 'learn-wisdom.txt')

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
      rename_file('wisdom.txt','learn-wisdom.txt')

      @wisdom_file.refresh
      expect(last_response.status).to eq(403)
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'refuses to rename a non-existent file' do
      # we attempt to rename a file that does not exist
      token_header(:editor)
      rename_file('folly.txt','learn-folly.txt')

      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('File not found')

      # the actual file is untouched
      @wisdom_file.refresh
      expect(@wisdom_file.file_name).to eq('wisdom.txt')
      expect(@wisdom_file).to be_has_data
    end

    it 'refuses to rename over an existing file' do
      learn_wisdom_file = create_file('athena', 'learn-wisdom.txt', WISDOM*2)
      stub_file('learn-wisdom.txt', WISDOM*2, :athena)

      token_header(:editor)
      rename_file('wisdom.txt','learn-wisdom.txt')

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
      expect(::File.read(@wisdom_file.location)).to eq(WISDOM)
      expect(::File.read(learn_wisdom_file.location)).to eq(WISDOM*2)
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
      stub_folder('contents', 'athena')

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
      stub_folder('contents', 'athena')

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

  context '#compute_hash!' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stub_file('wisdom.txt', WISDOM, :athena)
    end

    it 'computes the md5 sum of the file' do
      @wisdom_file.compute_hash!

      @wisdom_file.refresh
      expect(@wisdom_file.file_hash).to eq(Digest::MD5.hexdigest(WISDOM))
    end
  end

  context '#backup!' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stub_file('wisdom.txt', WISDOM, :athena)
    end

    def glacier_stub(method, path, response={})
      stub_request(method, "https://glacier.us-east-1.amazonaws.com/-/vaults#{path}")
        .to_return({
          status: 200,
          headers: {
          'Content-Type' => 'application/json'
          }
        }.merge(response))
    end

    it 'backs up the file to AWS Glacier' do
      vault_json = {
        "CreationDate":"2018-08-27T20:56:30.058Z",
        "LastInventoryDate":nil,
        "NumberOfArchives":0,
        "SizeInBytes":0,
        "VaultARN":"arn:aws:glacier:us-east-1:1999:vaults/metis-test-athena",
        "VaultName":"metis-test-athena"
      }.to_json

      glacier_stub( :get, '',
        body: '{"Marker":null,"VaultList":[]}',
      )
      glacier_stub(:put, '/metis-test-athena', status: 201, body: '{}')
      glacier_stub(:get, '/metis-test-athena',
        body: vault_json,
      )
      glacier_stub(:post, '/metis-test-athena',
        body: vault_json,
      )
      glacier_stub(:post, '/metis-test-athena/multipart-uploads',
        status: 201,
        headers: {
          'X-Amz-Multipart-Upload-Id': 'uploadid',
          'Content-Type': 'application/json'
        },
        body: "{}"
      )
      glacier_stub(:put, '/metis-test-athena/multipart-uploads/uploadid',
        status: 204
      )
      glacier_stub(:post, '/metis-test-athena/multipart-uploads/uploadid',
        status: 201,
        headers: {
          'X-Amz-Archive-Id': 'archive_id',
          'Content-Type': 'application/json',
        }
      )

      @wisdom_file.backup!

      @wisdom_file.refresh
      expect(@wisdom_file.archive_id).to eq('archive_id')
    end
  end
end
