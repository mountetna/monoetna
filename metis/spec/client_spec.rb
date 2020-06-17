load 'bin/metis_client'

def expect_output(path, *args)
  expect {
    begin
      MetisShell.new(path, *args).run
    rescue SystemExit
    end
  }.to output(yield).to_stdout
end

describe MetisShell do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    MetisConfig.instance.config(
      metis_uid: SecureRandom.hex,
      metis_uid_name: 'METIS_TEST_UID',
      metis_host: 'metis.test'
    )
    ENV['TOKEN'] = Base64.strict_encode64(
      { email: 'metis@olympus.org', perm: 'a:athena' }.to_json
    )
    stub_request(:any, %r!^https://metis.test/!).to_rack(app)
  end

  describe MetisShell::Ls do
    it 'lists root buckets' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      expect_output("metis://athena", "ls") { %r!armor/!}
    end

    it 'lists files and folders' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket)
      expect_output("metis://athena/armor", "ls") { %r!helmet.jpg.*helmet/!m }
    end

    it 'lists from a directory' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      expect_output("metis://athena/armor/helmet", "ls") { "helmet.jpg\n" }
    end

    it 'lists a parent directory' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      expect_output("metis://athena/armor/helmet", "ls", "..") { "helmet/\n" }
    end

    it 'lists the current directory' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      expect_output("metis://athena/armor/helmet", "ls", ".") { "helmet.jpg\n" }
    end

    it 'lists a sub-directory' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      expect_output("metis://athena/armor/", "ls", "helmet") { "helmet.jpg\n" }
    end
    it 'lists the project root directory' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      expect_output("metis://athena/armor/helmet", "ls", "/") { "armor/\n" }
    end
  end

  describe MetisShell::Put do
    xit 'puts a file into a bucket' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_path = stubs.create_data('stubs', 'helmet.txt', HELMET)
      #expect_output("metis://athena/armor", "put", helmet_path) { '' }
      #expect(Metis::File.first.file_name).to eq('helmet.txt')
    end
  end

  describe MetisShell::Get do
    xit 'downloads a folder from a bucket' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      # it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg/ }
      expect(::File.read("spec/data/helmet/helmet.jpg")).to eq(HELMET)

      # cleanup
      FileUtils.rm_rf('spec/data/helmet')
    end

    xit 'does not download data twice' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      # first time it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg/ }

      # second time it shows no progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { "" }

      # cleanup
      FileUtils.rm_rf('spec/data/helmet')
    end
  end
end
