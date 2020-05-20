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

    it 'returns the source bucket name' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.source_bucket).to eq('files')
    end

    it 'returns the dest bucket name' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.dest_bucket).to eq('files')
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

    it 'generates a full Metis path without folder_path' do
        path = Metis::Revision.path_from_parts('athena', 'files', 'learn-wisdom.txt')
        expect(path).to eq('metis://athena/files/learn-wisdom.txt')
    end

    it 'generates a full Metis path with folder_path' do
        path = Metis::Revision.path_from_parts('athena', 'files', 'helmet.jpg', 'blueprints/helmet')
        expect(path).to eq('metis://athena/files/blueprints/helmet/helmet.jpg')
    end
end