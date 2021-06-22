describe Metis::FileRenameRevision do

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

    @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
    stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  it 'creates a Metis::PathWithObjects as the dest parameter' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/wisdom.txt',
          user: @user
      })
      expect(revision.dest.instance_of? Metis::PathWithObjects).to eq(true)
  end

  it 'adds error message if user cannot access the dest bucket' do
      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'viewer', owner: 'metis' )

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/sundry/wisdom.txt',
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
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/wisdom2.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors).to eq([])
  end

  it 'adds error message if the dest path is invalid' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: "metis://athena/files/learn\nwisdom.txt",
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
          "Invalid path: \"metis://athena/files/learn\nwisdom.txt\""
      )

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
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
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/learn-wisdom.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors).to eq([])

      blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/blueprints/build-helmet.jpg',
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
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/helmet.jpg',
          dest: 'metis://athena/files/wisdom.txt',
          user: @user
      })
      expect(revision.bucket_names).
          to eq(['files', 'files'])

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/helmet.jpg',
          dest: nil,
          user: @user
      })
      expect(revision.bucket_names).
          to eq(['files'])
  end

  it 'adds error message if dest bucket is read-only' do
      contents_folder = create_folder('athena', 'contents', read_only: true)
      stubs.create_folder('athena', 'files', 'contents')

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/contents/wisdom.txt',
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

  it 'adds error message if source file exists and is read-only' do
      @wisdom2_file = create_file('athena', 'wisdom2.txt', WISDOM*2, read_only: true)
      stubs.create_file('athena', 'files', 'wisdom2.txt', WISDOM*2)

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom2.txt',
          dest: 'metis://athena/files/wisdom3.txt',
          user: @user
      })

      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
          "File \"metis://athena/files/wisdom2.txt\" is read-only"
      )
  end

  it 'adds error message if dest file exists' do
    @wisdom2_file = create_file('athena', 'wisdom2.txt', WISDOM*2)
    stubs.create_file('athena', 'files', 'wisdom2.txt', WISDOM*2)

    revision = Metis::FileRenameRevision.new({
        source: 'metis://athena/files/wisdom.txt',
        dest: 'metis://athena/files/wisdom2.txt',
        user: @user
    })

    revision.source.bucket = default_bucket('athena')
    revision.dest.bucket = default_bucket('athena')
    expect(revision.errors).to eq(nil)
    revision.validate
    expect(revision.errors.length).to eq(1)
    expect(revision.errors[0]).to eq(
        "File \"metis://athena/files/wisdom2.txt\" already exists"
    )
end

  it 'reports multiple errors from validation' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/helmet.jpg',
          dest: nil,
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate

      expect(revision.errors.length).to eq(2)

      expect(revision.errors[0]).to eq(
          "File \"metis://athena/files/helmet.jpg\" not found"
      )
      expect(revision.errors[1]).to eq(
          "Invalid path: \"\""
      )
  end

  it 'adds error message if trying to copy over an existing folder' do
      wisdom_folder = create_folder('athena', 'wisdom.txt')
      stubs.create_folder('athena', 'files', 'wisdom.txt')

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/wisdom.txt',
          user: @user
      })

      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
          "Cannot copy over existing folder: \"metis://athena/files/wisdom.txt\""
      )
  end

  it 'executes the revision' do
      expect(Metis::File.count).to eq(1)
      expect(Metis::File.first.author).to eq('metis|Metis')
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/learn-wisdom.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      revision.validate
      learn_wisdom = revision.revise!
      expect(Metis::File.count).to eq(1)
      expect(Metis::File.first.author).to eq(Metis::File.author(@user))
      expect(learn_wisdom.data_block).to eq(@wisdom_file.data_block)
  end

  it 'executes the revision with a new bucket' do
    new_bucket_name = 'new_bucket'
    stubs.create_bucket(new_bucket_name, 'files')
    new_bucket = create( :bucket, project_name: new_bucket_name, name: 'files', owner: 'metis', access: 'viewer')
    
    expect(Metis::File.count).to eq(1)
    expect(Metis::File.first.author).to eq('metis|Metis')
    revision = Metis::FileRenameRevision.new({
        source: 'metis://athena/files/wisdom.txt',
        dest: "metis://#{new_bucket_name}/files/learn-wisdom.txt",
        user: @user
    })
    revision.source.bucket = default_bucket('athena')
    revision.dest.bucket = new_bucket
    revision.validate
    learn_wisdom = revision.revise!
    expect(Metis::File.count).to eq(1)
    expect(Metis::File.first.author).to eq(Metis::File.author(@user))
    expect(learn_wisdom.data_block).to eq(@wisdom_file.data_block)
end

  it 'throws exception if you try to revise without setting user' do
      expect(Metis::File.count).to eq(1)
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/learn-wisdom.txt'
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      revision.validate
      expect {
          revision.revise!
      }.to raise_error(StandardError)

      expect(Metis::File.count).to eq(1)
  end

  it 'adds error message if source bucket is invalid' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/war/helmet.jpg',
          dest: 'metis://athena/files/helmet2.jpg',
          user: @user
      })
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
        "Invalid bucket \"war\" in project \"athena\". Check the bucket name and your permissions."
      )
  end

  it 'adds error message if source file does not exist' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/learn-wisdom.txt',
          dest: 'metis://athena/files/wisdom2.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
          "File \"metis://athena/files/learn-wisdom.txt\" not found"
      )
  end

  it 'adds error message if the source path is invalid' do
      revision = Metis::FileRenameRevision.new({
          source: "metis://athena/files/build\nhelmet.jpg",
          dest: "metis://athena/magma/learn-wisdom.txt",
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
          "Invalid path: \"metis://athena/files/build\nhelmet.jpg\""
      )

      revision = Metis::FileRenameRevision.new({
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
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/learn-wisdom.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors).to eq([])
  end

  it 'adds error message if user cannot access the source bucket' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/magma/wisdom.txt',
          dest: 'metis://athena/files/wisdom2.txt',
          user: @user
      })
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors.length).to eq(1)
      expect(revision.errors[0]).to eq(
        "Invalid bucket \"magma\" in project \"athena\". Check the bucket name and your permissions."
      )
  end

  it 'does not add error message if user can access the source bucket' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/magma/wisdom2.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      revision.validate
      expect(revision.errors).to eq([])
  end

  it 'will not execute the revision if not valid' do
      expect(Metis::File.count).to eq(1)
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
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

      expect(Metis::File.count).to eq(1)

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/wisdom.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      expect {
          revision.revise!
      }.to raise_error(StandardError)

      expect(Metis::File.count).to eq(1)
  end

  it 'is invalid if validate not run' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/wisdom.txt',
          dest: 'metis://athena/files/wisdom.txt',
          user: @user
      })
      revision.source.bucket = default_bucket('athena')
      revision.dest.bucket = default_bucket('athena')
      expect(revision.errors).to eq(nil)
      expect(revision.valid?).to eq(false)
  end

  it 'returns the dest and source paths in an array' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/helmet.jpg',
          dest: 'metis://athena/magma/wisdom.txt'
      })
      expect(revision.mpaths.length).to eq(2)
      expect(revision.mpaths[0].path).
          to eq('metis://athena/files/helmet.jpg')
      expect(revision.mpaths[1].path).
          to eq('metis://athena/magma/wisdom.txt')

      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/files/helmet.jpg',
          dest: nil
      })
      expect(revision.mpaths.length).to eq(1)
      expect(revision.mpaths[0].path).
          to eq('metis://athena/files/helmet.jpg')
  end

  it 'returns a hash representation' do
      revision = Metis::FileRenameRevision.new({
          source: 'metis://athena/magma/wisdom.txt',
          dest: 'metis://athena/files/wisdom2.txt',
          user: @user
      })
      revision.dest.bucket = default_bucket('athena')

      expect(revision.to_hash).to eq({
          source: 'metis://athena/magma/wisdom.txt',
          dest: 'metis://athena/files/wisdom2.txt',
          errors: nil
      })

      revision.validate

      expect(revision.to_hash).to eq({
          source: 'metis://athena/magma/wisdom.txt',
          dest: 'metis://athena/files/wisdom2.txt',
          errors: ["Invalid bucket \"magma\" in project \"athena\". Check the bucket name and your permissions."]
      })
  end
end
