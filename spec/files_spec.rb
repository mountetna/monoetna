describe FilesController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    clear_stubs
  end

  context '#index' do
    it 'should return a list of files' do
      # our files
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stub_file('wisdom.txt', WISDOM, :athena)

      helmet = 'x'*20
      helmet_file = create_file('athena', 'helmet.jpg', helmet)
      stub_file('helmet.jpg', helmet, :athena)

      # we post a cancel request with our hmac url
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      get('/athena/files')

      expect(last_response.status).to eq(200)
      json = json_body(last_response.body)
      expect(json).to eq(files: [
        {file_name: "wisdom.txt", project_name: "athena", original_name: "wisdom.txt", size: 66, file_hash: Digest::MD5.hexdigest(WISDOM)},
        {file_name: "helmet.jpg", project_name: "athena", original_name: "helmet.jpg", size: 20, file_hash: Digest::MD5.hexdigest(helmet)}
      ])
    end
  end
end
