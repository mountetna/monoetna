describe Metis::Path do
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

      expect(stubs.contents(:athena)).to be_empty
    end

    it 'creates correctly' do
        path = Metis::Path.new('metis://athena/files/helmet.jpg')
    end

    it 'returns the bucket name' do
        path = Metis::Path.new('metis://athena/files/helmet.jpg')
        expect(path.bucket.name).to eq('files')
    end

    it 'returns a full Metis path from parts' do
        path = Metis::Path.path_from_parts('athena', 'files', 'learn-wisdom.txt')
        expect(path).to eq('metis://athena/files/learn-wisdom.txt')

        path = Metis::Path.path_from_parts('athena', 'files', 'blueprints/helmet/helmet.jpg')
        expect(path).to eq('metis://athena/files/blueprints/helmet/helmet.jpg')
    end

    it 'correctly returns the file path' do
        path = Metis::Path.new(
            'metis://athena/files/blueprints/helmet/helmet.jpg'
        )
        expect(path.file_path).to eq('blueprints/helmet/helmet.jpg')
    end

    it 'returns false for invalid path' do
        path = Metis::Path.new(
                "metis://athena/files/learn\nwisdom.txt")
        expect(path.valid?).to eq(false)

        path = Metis::Path.new(nil)
        expect(path.valid?).to eq(false)
    end

    it 'returns true for valid path' do
        path = Metis::Path.new(
                "metis://athena/files/learn-wisdom.txt")
        expect(path.valid?).to eq(true)
    end

    it 'returns the path\'s file' do
        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        path = Metis::Path.new('metis://athena/files/wisdom.txt')
        expect(path.file).to eq(@wisdom_file)
    end

    it 'returns the path\'s folder' do
        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        path = Metis::Path.new('metis://athena/files/wisdom.txt')
        expect(path.folder).to eq(@wisdom_file.folder)
    end

    it 'returns the path\'s folder_path' do
        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        path = Metis::Path.new('metis://athena/files/wisdom.txt')
        expect(path.folder_path).to eq(nil)

        path = Metis::Path.new('metis://athena/files/blueprints/helmet/helmet.jpg')
        expect(path.folder_path).to eq('blueprints/helmet')
    end

    it 'returns the path\'s file_name' do
        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        path = Metis::Path.new('metis://athena/files/wisdom.txt')
        expect(path.file_name).to eq('wisdom.txt')

        path = Metis::Path.new('metis://athena/files/blueprints/helmet/helmet.jpg')
        expect(path.file_name).to eq('helmet.jpg')
    end

    it 'returns the path\'s bucket' do
        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        path = Metis::Path.new('metis://athena/files/wisdom.txt')
        expect(path.bucket).to eq(@wisdom_file.bucket)
    end

    it 'returns the path\'s bucket name' do
        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        path = Metis::Path.new('metis://athena/files/wisdom.txt')
        expect(path.bucket_name).to eq(@wisdom_file.bucket.name)
    end
end