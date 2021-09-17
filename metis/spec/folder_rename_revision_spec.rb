describe Metis::FolderRenameRevision do

    def app
      OUTER_APP
    end

    before(:each) do
      default_bucket('athena')

      @user = Etna::User.new({
        name: 'Athena',
        email: 'athena@olympus.org',
        perm: 'a:athena'
      })

      @wisdom_folder = create_folder('athena', 'wisdom')
      stubs.create_folder('athena', 'files', 'wisdom')
    end

    after(:each) do
      stubs.clear

      expect(stubs.contents(:athena)).to be_empty
    end

    it 'creates a Metis::PathWithObjects as the dest parameter' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/wisdom',
            user: @user
        })
        expect(revision.dest.instance_of? Metis::PathWithObjects).to eq(true)
    end

    it 'adds error message if user cannot access the dest bucket' do
        sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/sundry/wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid bucket \"sundry\" in project \"athena\". Check the bucket name and your permissions."
        )
    end

    it 'does not add error message if user can access the dest bucket' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/wisdom_2',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors).to eq([])
    end

    it 'adds error message if the dest path is invalid' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: "metis://athena/files/learn\nwisdom",
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid path: \"metis://athena/files/learn\nwisdom\""
        )

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: nil,
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid path: \"\""
        )
    end

    it 'does not add error message if the dest path is valid' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/learn-wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors).to eq([])

        blueprints_folder = create_folder('athena', 'blueprints')
        stubs.create_folder('athena', 'files', 'blueprints')

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/blueprints/drawing_wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        revision.dest.folder = blueprints_folder
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors).to eq([])
    end

    it 'returns the source and dest bucket_names in an array' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/helmet',
            dest: 'metis://athena/files/wisdom',
            user: @user
        })
        expect(revision.bucket_names).
            to eq(['files', 'files'])

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/helmet',
            dest: nil,
            user: @user
        })
        expect(revision.bucket_names).
            to eq(['files'])
    end

    it 'adds error message if dest bucket is read-only' do
        contents_folder = create_folder('athena', 'contents', read_only: true)
        stubs.create_folder('athena', 'files', 'contents')

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/contents/wisdom',
            user: @user
        })

        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        revision.dest.folder = contents_folder
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Folder \"contents\" is read-only"
        )
    end

    it 'adds error message if dest file exists' do
        @wisdom2_file = create_file('athena', 'wisdom2.txt', WISDOM*2)
        stubs.create_file('athena', 'files', 'wisdom2.txt', WISDOM*2)

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/wisdom2.txt',
            user: @user
        })

        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Cannot overwrite existing file: \"metis://athena/files/wisdom2.txt\""
        )
    end

    it 'reports multiple errors from validation' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/helmet',
            dest: nil,
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate

        expect(revision.errors.length).to eq(2)

        expect(revision.errors[0]).to eq(
            "Folder not found: \"metis://athena/files/helmet\""
        )
        expect(revision.errors[1]).to eq(
            "Invalid path: \"\""
        )
    end

    it 'adds error message if trying to copy over an existing folder' do
        wisdom_folder = create_folder('athena', 'wisdom_jr')
        stubs.create_folder('athena', 'files', 'wisdom_jr')

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/wisdom_jr',
            user: @user
        })

        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Cannot copy over existing folder: \"metis://athena/files/wisdom_jr\""
        )
    end

    it 'executes the revision' do
        expect(Metis::Folder.count).to eq(1)
        expect(Metis::Folder.first.author).to eq('metis|Metis')
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/learn-wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        revision.validate
        learn_wisdom = revision.revise!
        expect(Metis::Folder.count).to eq(1)
        expect(Metis::Folder.first.author).to eq(Metis::File.author(@user))
        expect(learn_wisdom.id).to eq(@wisdom_folder.id)
    end

    it 'executes the revision to a new bucket' do
        new_bucket_name = 'new_bucket'
        stubs.create_bucket(new_bucket_name, 'files')
        new_bucket = create( :bucket, project_name: new_bucket_name, name: 'files', owner: 'metis', access: 'viewer')
        
        expect(Metis::Folder.count).to eq(1)
        expect(Metis::Folder.first.author).to eq('metis|Metis')
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: "metis://#{new_bucket_name}/files/learn-wisdom",
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = new_bucket
        revision.validate
        learn_wisdom = revision.revise!
        expect(Metis::Folder.count).to eq(1)
        expect(Metis::Folder.first.author).to eq(Metis::File.author(@user))
        expect(learn_wisdom.id).to eq(@wisdom_folder.id)
    end

    it 'throws exception if you try to revise without setting user' do
        expect(Metis::Folder.count).to eq(1)
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/learn-wisdom'
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        revision.validate
        expect {
            revision.revise!
        }.to raise_error(StandardError)

        expect(Metis::Folder.count).to eq(1)
    end

    it 'adds error message if source bucket is invalid' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/war/helmet',
            dest: 'metis://athena/files/wisdom',
            user: @user
        })
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(2)
        expect(revision.errors).to eq([
            "Invalid bucket \"war\" in project \"athena\". Check the bucket name and your permissions.",
            "Cannot copy over existing folder: \"metis://athena/files/wisdom\""
        ])
    end

    it 'adds error message if source folder does not exist' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom_jr',
            dest: 'metis://athena/files/wisdom_sr',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Folder not found: \"metis://athena/files/wisdom_jr\""
        )
    end

    it 'adds error message if the source path is invalid' do
        revision = Metis::FolderRenameRevision.new({
            source: "metis://athena/files/build\nhelmet",
            dest: "metis://athena/magma/instructables",
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect(revision.errors[0]).to eq(
            "Invalid path: \"metis://athena/files/build\nhelmet\""
        )

        revision = Metis::FolderRenameRevision.new({
            source: nil,
            dest: nil,
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(2)
        expect(revision.errors[0]).to eq(
            "Invalid path: \"\""
        )
        expect(revision.errors[1]).to eq(
            "Invalid path: \"\""
        )
    end

    it 'does not add error message if the source path is valid' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/backup_wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors).to eq([])
    end

    it 'adds error message if user cannot access the source bucket' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/magma/wisdom',
            dest: 'metis://athena/files/wisdom',
            user: @user
        })
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(2)
        expect(revision.errors[0]).to eq(
            "Invalid bucket \"magma\" in project \"athena\". Check the bucket name and your permissions."
        )
    end

    it 'does not add error message if user can access the source bucket' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/wisdom_jr',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors).to eq([])
    end

    it 'will not execute the revision if not valid' do
        expect(Metis::Folder.count).to eq(1)
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: nil,
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        revision.validate
        expect(revision.errors.length).to eq(1)
        expect {
            revision.revise!
        }.to raise_error(StandardError)

        expect(Metis::Folder.count).to eq(1)

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/files/backup_wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        expect {
            revision.revise!
        }.to raise_error(StandardError)

        expect(Metis::Folder.count).to eq(1)
    end

    it 'is invalid if validate not run' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/wisdom',
            dest: 'metis://athena/backup_files/wisdom',
            user: @user
        })
        revision.source.bucket = default_bucket('athena')
        revision.dest.bucket = default_bucket('athena')
        expect(revision.errors).to eq(nil)
        expect(revision.valid?).to eq(false)
    end

    it 'returns the dest and source paths in an array' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/helmet',
            dest: 'metis://athena/magma/wisdom'
        })
        expect(revision.mpaths.length).to eq(2)
        expect(revision.mpaths[0].path).
            to eq('metis://athena/files/helmet')
        expect(revision.mpaths[1].path).
            to eq('metis://athena/magma/wisdom')

        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/files/helmet',
            dest: nil
        })
        expect(revision.mpaths.length).to eq(1)
        expect(revision.mpaths[0].path).
            to eq('metis://athena/files/helmet')
    end

    it 'returns a hash representation' do
        revision = Metis::FolderRenameRevision.new({
            source: 'metis://athena/magma/wisdom',
            dest: 'metis://athena/files/wisdom',
            user: @user
        })
        revision.dest.bucket = default_bucket('athena')

        expect(revision.to_hash).to eq({
            source: 'metis://athena/magma/wisdom',
            dest: 'metis://athena/files/wisdom',
            errors: nil
        })

        revision.validate

        expect(revision.to_hash).to eq({
            source: 'metis://athena/magma/wisdom',
            dest: 'metis://athena/files/wisdom',
            errors: [
                "Invalid bucket \"magma\" in project \"athena\". Check the bucket name and your permissions.",
                "Cannot copy over existing folder: \"metis://athena/files/wisdom\""],
        })
    end
end
