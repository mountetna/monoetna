describe Metis::CopyRevision do
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

    it 'creates a CopyRevision from parts' do
        revision = Metis::CopyRevision.create_from_parts({
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
        expect(revision.instance_of? Metis::CopyRevision).to be_truthy
    end

    it 'creates a Metis::Path as the dest parameter' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.dest.instance_of? Metis::Path).to eq(true)
    end

    it 'returns false if user cannot access the dest bucket' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.valid?('dest_bucket_access', ['files'])).
          to eq(false)
    end

    it 'returns true if user can access the dest bucket' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.valid?('dest_bucket_access', ['magma'])).
            to eq(true)
    end

    it 'returns false if the dest path is invalid' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: "metis://athena/magma/learn\nwisdom.txt"
        })
        expect(revision.valid?('dest_path')).
          to eq(false)

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: nil
        })
        expect(revision.valid?('dest_path')).
          to eq(false)
    end

    it 'returns true if the dest path is valid' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.valid?('dest_path')).
            to eq(true)
    end

    it 'returns the source and dest bucket_names in an array' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.bucket_names).
            to eq(['files', 'magma'])
    end
end