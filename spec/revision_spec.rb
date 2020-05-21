describe Metis::Revision do
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

    it 'throws exception if no source provided' do
        expect {
            Metis::Revision.new({
                dest: 'metis://athena/files/wisdom.txt'
            })
        }.to raise_error(Etna::BadRequest)
    end

    it 'throws exception if no dest provided' do
        expect {
            Metis::Revision.new({
                source: 'metis://athena/files/wisdom.txt'
            })
        }.to raise_error(Etna::BadRequest)
    end

    it 'creates correctly' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
    end

    it 'creates correctly with nil dest' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: nil
        })
    end

    it 'returns the source bucket name' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.source_bucket_name).to eq('files')
    end

    it 'returns the dest bucket name' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.dest_bucket_name).to eq('files')
    end

    it 'throws exception if user cannot access the source bucket' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect {
            revision.validate_access_to_buckets(['magma'])
        }.to raise_error(Etna::Forbidden)
    end

    it 'throws exception if user cannot access the destination bucket' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect {
            revision.validate_access_to_buckets(['files'])
        }.to raise_error(Etna::Forbidden)
    end

    it 'does nothing if user can access the source and dest buckets' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate_access_to_buckets(['magma', 'files'])
    end

    it 'generates a full Metis path from parts' do
        path = Metis::Revision.path_from_parts('athena', 'files', 'learn-wisdom.txt')
        expect(path).to eq('metis://athena/files/learn-wisdom.txt')

        path = Metis::Revision.path_from_parts('athena', 'files', 'blueprints/helmet/helmet.jpg')
        expect(path).to eq('metis://athena/files/blueprints/helmet/helmet.jpg')
    end

    it 'returns the bucket for a given path' do
        expect(Metis::Revision.extract_bucket_name_from_path(
            'metis://athena/files/blueprints/helmet/helmet.jpg'
        )).to eq('files')
    end

    it 'returns the file_path for a given path' do
        expect(Metis::Revision.extract_file_path_from_path(
            'metis://athena/files/blueprints/helmet/helmet.jpg'
        )).to eq('blueprints/helmet/helmet.jpg')

        expect(Metis::Revision.extract_file_path_from_path(
            'metis://athena/files/wisdom.txt'
        )).to eq('wisdom.txt')
    end

    it 'correctly returns the source file path' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.source_file_path).to eq('blueprints/helmet/helmet.jpg')
    end

    it 'correctly returns the dest file path' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/blueprints/helmet/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.dest_file_path).to eq('wisdom.txt')
    end

    it 'throws exception for invalid source path' do
        expect {
            Metis::Revision.new({
                source: "metis://athena/files/learn\nwisdom.txt",
                dest: nil
            }).to raise_error(Etna::BadRequest)
        }
    end
end