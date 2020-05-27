describe Metis::CopyRevision do
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
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.dest.instance_of? Metis::Path).to eq(true)
    end

    it 'adds error message if user cannot access the dest bucket' do
        sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/sundry/wisdom.txt'
        })
        revision.validate(['files'])

        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Forbidden: no access to dest bucket sundry"
        )
    end

    it 'does not add error message if user can access the dest bucket' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/wisdom.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(0)
    end

    it 'adds error message if the dest path is invalid' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: "metis://athena/files/learn\nwisdom.txt"
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid dest path: metis://athena/files/learn\nwisdom.txt"
        )

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: nil
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid dest path: "
        )
    end

    it 'does not add error message if the dest path is valid' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/learn-wisdom.txt'
        })
        revision.validate(['files'])

        expect(revision.errors.length).to eq(0)

        blueprints_folder = create_folder('athena', 'blueprints')
        stubs.create_folder('athena', 'files', 'blueprints')

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/blueprints/build-helmet.jpg'
        })
        revision.validate(['files'])

        expect(revision.errors.length).to eq(0)
    end

    it 'returns the source and dest bucket_names in an array' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.bucket_names).
            to eq(['files', 'files'])

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: nil
        })
        expect(revision.bucket_names).
            to eq(['files'])
    end

    it 'adds error message if dest bucket is read-only' do
        contents_folder = create_folder('athena', 'contents', read_only: true)
        stubs.create_folder('athena', 'files', 'contents')

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/contents/wisdom.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Dest folder contents is read-only"
        )
    end

    it 'adds error message if dest file exists and is read-only' do
        @wisdom2_file = create_file('athena', 'wisdom2.txt', WISDOM*2, read_only: true)
        stubs.create_file('athena', 'files', 'wisdom2.txt', WISDOM*2)

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/wisdom2.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Dest file wisdom2.txt is read-only"
        )
    end

    it 'reports multiple errors from validation' do
        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: nil
        })
        revision.validate(['magma'])
        expect(revision.errors.length).to eq(3)

        expect(revision.errors[0]).to eq(
            "File metis://athena/files/helmet.jpg not found"
        )
        expect(revision.errors[1]).to eq(
            "Forbidden: no access to source bucket files"
        )
        expect(revision.errors[2]).to eq(
            "Invalid dest path: "
        )
    end

    it 'adds error message if trying to copy over an existing bucket' do
        wisdom_folder = create_folder('athena', 'wisdom.txt')
        stubs.create_folder('athena', 'files', 'wisdom.txt')

        revision = Metis::CopyRevision.new({
            source: 'metis://athena/files/wisdom.txt',
            dest: 'metis://athena/files/wisdom.txt'
        })
        revision.validate(['files'])
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Cannot copy over existing folder metis://athena/files/wisdom.txt"
        )
    end
end