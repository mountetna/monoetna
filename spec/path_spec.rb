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
        expect(path.bucket_name).to eq('files')
    end

    it 'returns a full Metis path from parts' do
        path = Metis::Path.path_from_parts('athena', 'files', 'learn-wisdom.txt')
        expect(path).to eq('metis://athena/files/learn-wisdom.txt')

        path = Metis::Path.path_from_parts('athena', 'files', 'blueprints/helmet/helmet.jpg')
        expect(path).to eq('metis://athena/files/blueprints/helmet/helmet.jpg')
    end

    it 'returns the bucket name for a given path' do
        expect(Metis::Path.extract_bucket_name_from_path(
            'metis://athena/files/blueprints/helmet/helmet.jpg'
        )).to eq('files')
    end

    it 'returns the file_path for a given path' do
        expect(Metis::Path.extract_file_path_from_path(
            'metis://athena/files/blueprints/helmet/helmet.jpg'
        )).to eq('blueprints/helmet/helmet.jpg')

        expect(Metis::Path.extract_file_path_from_path(
            'metis://athena/files/wisdom.txt'
        )).to eq('wisdom.txt')
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
end