describe Metis::Revision do
    def app
      OUTER_APP
    end

    before(:each) do
      default_bucket('athena')

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

    it 'creates a Metis::PathWithObjects as the source parameter' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/files/wisdom.txt'
        })
        expect(revision.source.instance_of? Metis::PathWithObjects).to eq(true)
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

    it 'returns the source path in an array' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.mpaths.length).to eq(1)
        expect(revision.mpaths[0].path).
            to eq('metis://athena/files/helmet.jpg')

        revision = Metis::Revision.new({
            source: nil,
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.mpaths).
            to eq([])
    end

    it 'sets the bucket on a Metis Path with Objects' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.source.bucket).to eq(nil)
        revision.set_bucket(revision.source, [default_bucket('athena')])
        expect(revision.source.bucket).to eq(default_bucket('athena'))
    end

    it 'cannot set bucket if invalid Metis Path' do
        revision = Metis::Revision.new({
            source: 'metis://athena/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.source.bucket).to eq(nil)
        revision.set_bucket(revision.source, [default_bucket('athena')])
        expect(revision.source.bucket).to eq(nil)
    end

    it 'does not set the bucket on a Metis Path with Objects if does not match' do
        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.dest.bucket).to eq(nil)
        revision.set_bucket(revision.dest, [default_bucket('athena')])
        expect(revision.dest.bucket).to eq(nil)
    end

    it 'sets the folder on a Metis Path with Objects' do
        blueprints_folder = create_folder('athena', 'blueprints')
        stubs.create_folder('athena', 'files', 'blueprints')

        revision = Metis::Revision.new({
            source: 'metis://athena/files/blueprints/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.set_bucket(revision.source, [default_bucket('athena')])

        expect(revision.source.folder).to eq(nil)
        revision.set_folder(revision.source, [blueprints_folder])
        expect(revision.source.folder).to eq(blueprints_folder)
    end

    it 'does not the folder on a Metis Path with Objects if does not match path' do
        blueprints_folder = create_folder('athena', 'blueprints')
        stubs.create_folder('athena', 'files', 'blueprints')

        revision = Metis::Revision.new({
            source: 'metis://athena/files/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.set_bucket(revision.source, [default_bucket('athena')])

        expect(revision.source.folder).to eq(nil)
        revision.set_folder(revision.source, [blueprints_folder])
        expect(revision.source.folder).to eq(nil)
    end

    it 'cannot set folder if invalid Metis Path' do
        blueprints_folder = create_folder('athena', 'blueprints')
        stubs.create_folder('athena', 'files', 'blueprints')

        revision = Metis::Revision.new({
            source: 'metis://athena/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        revision.set_bucket(revision.source, [default_bucket('athena')])
        expect(revision.source.folder).to eq(nil)
        revision.set_folder(revision.source, [blueprints_folder])
        expect(revision.source.folder).to eq(nil)
    end

    it 'cannot set folder if no bucket on Metis Path' do
        blueprints_folder = create_folder('athena', 'blueprints')
        stubs.create_folder('athena', 'files', 'blueprints')

        revision = Metis::Revision.new({
            source: 'metis://athena/helmet.jpg',
            dest: 'metis://athena/magma/wisdom.txt'
        })
        expect(revision.source.folder).to eq(nil)
        revision.set_folder(revision.source, [blueprints_folder])
        expect(revision.source.folder).to eq(nil)
    end
end