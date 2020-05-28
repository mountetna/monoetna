describe Metis::Revision do
    include Rack::Test::Methods

    def app
      OUTER_APP
    end

    before(:each) do
      default_bucket('athena')

      @metis_uid = Metis.instance.sign.uid

      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"

      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
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

    it 'throws exception if you try to directly call the abstract validate method' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect {
            revision.validate([])
        }.to raise_error(StandardError)
    end

    it 'adds error message if user cannot access the source bucket' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate(['magma'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Forbidden: no access to source bucket files"
        )
    end

    it 'does not add error message if user can access the source bucket' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(0)
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

    it 'adds error message if the source path is invalid' do
        revision = Metis::Revision.new({
            source: "metis://athena/files/build\nhelmet.jpg",
            dest: "metis://athena/magma/learn-wisdom.txt"
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid source path: metis://athena/files/build\nhelmet.jpg"
        )

        revision = Metis::Revision.new({
            source: nil,
            dest: nil
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid source path: "
        )
    end

    it 'does not add error message if the source path is valid' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(0)
    end

    it 'returns the source bucket_name in an array' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.bucket_names).
            to eq(['files'])

        revision = Metis::Revision.new({
            source: nil,
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.bucket_names).
            to eq([])
    end

    it 'adds error message if source bucket is invalid' do
        revision = Metis::Revision.new({
            source: 'metis://athena/war/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate(['war'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid source bucket: war"
        )
    end

    it 'adds error message if source file does not exist' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/learn-wisdom.txt',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "File metis://athena/files/learn-wisdom.txt not found"
        )
    end

    it 'reports multiple errors from validation' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.validate(['magma'])
        expect(revision.errors.length).to eq(2)

        expect(revision.errors[0]).to eq(
            "File metis://athena/files/helmet.jpg not found"
        )
        expect(revision.errors[1]).to eq(
            "Forbidden: no access to source bucket files"
        )
    end
end