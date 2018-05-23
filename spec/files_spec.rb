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

  context '#list' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stub_file('wisdom.txt', WISDOM, :athena)

      @blueprints_folder = create_folder('athena', 'blueprints')

      @helmet = 'x'*20
      @helmet_file = create_file('athena', 'blueprints/helmet.jpg', @helmet, folder: @blueprints_folder)
      stub_file('blueprints/helmet.jpg', @helmet, :athena)
    end

    it 'should return a list of files and folders for the current folder' do
      # our files
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      get('/athena/list/files/')

      expect(last_response.status).to eq(200)

      expect(json_body[:files][0]).to include(
        file_name: 'wisdom.txt',
        author: 'metis|Metis',
        project_name: 'athena',
        size: 66,
        file_hash: Digest::MD5.hexdigest(WISDOM),
        download_url: a_string_matching(%r{http.*athena/download})
      )
      expect(json_body[:files][1]).to include(
        file_name: 'blueprints',
        author: 'metis|Metis',
        project_name: 'athena',
        is_folder: true
      )
    end

    it 'should list files from a sub-folder' do
      # our files
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      get('/athena/list/files/blueprints')

      expect(last_response.status).to eq(200)

      expect(json_body[:files].first).to include(
        file_name: 'blueprints/helmet.jpg',
        author: 'metis|Metis',
        project_name: 'athena',
        size: 20,
        file_hash: Digest::MD5.hexdigest(@helmet),
        download_url: a_string_matching(%r{http.*athena/download})
      )
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
      expect(json_body[:error]).to eq('Invalid folder name')
    end

    it 'creates nested folders' do
      blueprints_folder = create_folder('athena', 'blueprints')
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
      expect(json_body[:error]).to eq('Invalid parent folder')
    end

    it 'refuses to create folders with non-existent parent folder' do
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('athena/create_folder/files/blueprints/Helmet Blueprints', {})

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid parent folder')
    end

    it 'sets a parent folder' do
      blueprints_folder = create_folder('athena', 'blueprints')
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('athena/create_folder/files/blueprints/Helmet Blueprints', {})

      file = Metis::File.last

      expect(last_response.status).to eq(200)
      expect(Metis::File.count).to eq(2)
      expect(file.file_name).to eq('blueprints/Helmet Blueprints')
      expect(file).to be_folder
      expect(file.folder).to eq(blueprints_folder)
    end
  end
end
