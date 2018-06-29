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

      @helmet_folder = create_folder('athena', 'helmet', folder: @blueprints_folder)

      @helmet = 'x'*20
      @helmet_file = create_file('athena', 'helmet.jpg', @helmet, folder: @helmet_folder)
      stub_file('blueprints/helmet/helmet.jpg', @helmet, :athena)
    end

    it 'should return a list of files and folders for the current folder' do
      # our files
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      get('/athena/list/files/')

      expect(last_response.status).to eq(200)

      expect(json_body[:files].first).to include(
        file_name: 'wisdom.txt',
        author: 'metis|Metis',
        project_name: 'athena',
        bucket_name: 'files',
        size: 66,
        file_hash: Digest::MD5.hexdigest(WISDOM),
        download_url: a_string_matching(%r{http.*athena/download})
      )
      expect(json_body[:folders].first).to include(
        folder_name: 'blueprints',
        author: 'metis|Metis',
        project_name: 'athena',
        bucket_name: 'files'
      )
    end

    it 'should list files from a sub-folder' do
      # our files
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      get('/athena/list/files/blueprints/helmet')

      expect(last_response.status).to eq(200)

      expect(json_body[:files].first).to include(
        file_name: 'helmet.jpg',
        author: 'metis|Metis',
        project_name: 'athena',
        size: @helmet.length,
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

      folder = Metis::Folder.first
      expect(folder).not_to be_nil
      expect(folder.folder_name).to eq('Helmet Blueprints')
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

      folder = Metis::Folder.last
      expect(folder).not_to be_nil
      expect(folder.folder_path).to eq(['blueprints', 'Helmet Blueprints'])
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

      folder = Metis::Folder.last

      expect(last_response.status).to eq(200)
      expect(Metis::Folder.count).to eq(2)
      expect(folder.folder_path).to eq([ 'blueprints', 'Helmet Blueprints'])
      expect(folder.folder).to eq(blueprints_folder)
    end
  end
end
