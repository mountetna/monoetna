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
  end

  describe MetisShell::Put do
    it 'puts a file into a bucket' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_path = stubs.create_data('stubs', 'helmet.txt', HELMET)
      #expect_output("metis://athena/armor", "put", helmet_path) { '' }
      #expect(Metis::File.first.file_name).to eq('helmet.txt')
    end
  end
end
