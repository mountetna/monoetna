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

    it 'creates a Metis::Path as the source parameter' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.source.instance_of? Metis::Path).to eq(true)
    end

    it 'returns false if user cannot access the source bucket' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.valid?('source_bucket_access', ['magma'])).
          to eq(false)
    end

    it 'returns true if user can access the source bucket' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.valid?('source_bucket_access', ['files'])).
          to eq(true)
    end

    it 'creates a Revision from parts' do
        revision = Metis::Revision.create_from_parts({
            source: {
                project_name: 'athena',
                bucket_name: 'files',
                file_path: 'blueprints/helmet/helmet.jpg'
            },
            dest: {
                project_name: 'athena',
                bucket_name: 'files',
                file_path: 'build-helmet.jpg'
            }
        })

        expect(revision.source.path).to eq('metis://athena/files/blueprints/helmet/helmet.jpg')
        expect(revision.dest.path).to eq('metis://athena/files/build-helmet.jpg')
    end

    it 'returns false if the source path is invalid' do
        revision = Metis::Revision.new({
            source: "metis://athena/files/build\nhelmet.jpg",
            dest: "metis://athena/magma/learn-wisdom.txt"
        })
        expect(revision.valid?('source_path')).
          to eq(false)

        revision = Metis::Revision.new({
            source: nil,
            dest: nil
        })
        expect(revision.valid?('source_path')).
          to eq(false)
    end

    it 'returns true if the source path is valid' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.valid?('source_path')).
            to eq(true)
    end

    it 'returns the source bucket_name in an array' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.bucket_names).
            to eq(['files'])
    end
end