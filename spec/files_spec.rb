describe FilesController do
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

      # isolate the urls since they are harder to match
      urls = json[:files].map{|f| f.delete(:download_url)}

      expect(json).to eq(files: [
        {file_name: "wisdom.txt", project_name: "athena", size: 66, file_hash: Digest::MD5.hexdigest(WISDOM)},
        {file_name: "helmet.jpg", project_name: "athena", size: 20, file_hash: Digest::MD5.hexdigest(helmet)}
      ])

      expect(urls).to all( match(%r{http.*athena/download}) )
    end
  end

  context '#create_folder' do
    it 'creates a folder with the given name' do
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('athena/create_folder/files/Helmet Blueprints', {})

      expect(last_response.status).to eq(200)

      folder = Metis::File.first
      expect(folder).not_to be_nil
      expect(folder.file_name).to eq('Helmet Blueprints')
      expect(folder).to be_folder
    end

    it 'refuses to create folders with invalid names' do
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post("athena/create_folder/files/Helmet\nBlueprints", {})

      expect(last_response.status).to eq(422)
      json = json_body(last_response.body)
      expect(json[:error]).to eq('Invalid folder name')
    end

    it 'creates nested folders' do
      blueprints_folder = create_file('athena', 'blueprints', nil, is_folder: true)
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('athena/create_folder/files/blueprints/Helmet Blueprints', {})

      expect(last_response.status).to eq(200)

      folder = Metis::File.last
      expect(folder).not_to be_nil
      expect(folder.file_name).to eq('blueprints/Helmet Blueprints')
      expect(folder).to be_folder
    end

    it 'refuses to set a file as parent' do
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)

      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('athena/create_folder/files/wisdom.txt/Helmet Blueprints', {})

      expect(last_response.status).to eq(422)
      json = json_body(last_response.body)
      expect(json[:error]).to eq('Invalid parent folder')
    end

    it 'refuses to create folders with no parent folder' do
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('athena/create_folder/files/blueprints/Helmet Blueprints', {})

      expect(last_response.status).to eq(422)
      json = json_body(last_response.body)
      expect(json[:error]).to eq('Invalid parent folder')
    end
  end
end
