load 'bin/metis_client'

def expect_output(path, *args)
  expect {
    begin
      MetisShell.new(path, *args).run
    rescue SystemExit
    end
  }.to output(yield).to_stdout
end

def expect_error(path, *args)
  expect {
    begin
      MetisShell.new(path, *args).run
    rescue SystemExit
    end
  }.to output(yield).to_stderr_from_any_process
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

    it 'lists files and folders in long format' do
      Timecop.freeze(DateTime.parse("2020-06-17T04:37"))
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket)
      stubs.create_file('athena', 'armor', 'helmet.jpg', HELMET)

      expect_output("metis://athena/armor", "ls", "-l") {
        "metis    Jun 17 04:37    helmet/\n"+
        "metis 13 Jun 17 04:37 helmet.jpg\n"
      }
      Timecop.return
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
    # WebMock does not handle multipart/form-data, so it
    #   strips out the body data on those POST requests.
    # So this is actually impossible to test except in
    #   a 'production' environment.
    xit 'puts a file into a bucket' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_path = stubs.create_data('stubs', 'helmet.txt', HELMET)
      #expect_output("metis://athena/armor", "put", helmet_path) { '' }
      #expect(Metis::File.first.file_name).to eq('helmet.txt')
    end

    it 'retries to upload a file into a bucket on a 422 response' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_path = stubs.create_data('stubs', 'helmet.txt', HELMET)

      expect_output("metis://athena/armor", "put", helmet_path, '.') { nil }

      # There is a call to authorize the upload
      expect(WebMock).to have_requested(:post, /https:\/\/metis.test\/authorize\/upload/).
        with(body: hash_including({
          "project_name": "athena",
          "bucket_name": "armor",
          "file_path": "helmet.txt"
        }))

      # There is a non-reset call to start the upload
      expect(WebMock).to have_requested(:post, /https:\/\/metis.test\/athena\/upload\/armor\/helmet.txt/).
        with(query: hash_including({
          "X-Etna-Id": "metis"
        })).
        with(headers: {
          "Content-Type": "application/json"
        }).
        with { |req| (req.body.include? 'start') && !(req.body.include? 'reset')}

      # There are two attempts to upload blobs. One for the original attempt,
      #   and one for the retry.
      expect(WebMock).to have_requested(:post, /https:\/\/metis.test\/athena\/upload\/armor\/helmet.txt/).
        with(query: hash_including({
          "X-Etna-Id": "metis"
        })).
        with(headers: {
          "Content-Type": "multipart/form-data"
        }).twice

      # There is a call to reset the upload
      expect(WebMock).to have_requested(:post, /https:\/\/metis.test\/athena\/upload\/armor\/helmet.txt/).
        with(query: hash_including({
          "X-Etna-Id": "metis"
        })).
        with(headers: {
          "Content-Type": "application/json"
        }).
        with { |req| (req.body.include? 'start') && (req.body.include? 'reset')}
    end
  end

  describe MetisShell::Get do
    it 'downloads a folder from a bucket' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      # it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg/ }

      # This ideally would equal HELMET, but because we don't use Apache X-Sendfile
      #   while testing, no bits are sent by Metis and the local file has no data.
      # expect(::File.read("spec/data/helmet/helmet.jpg")).to eq(HELMET)

      expect(::File.read("spec/data/helmet/helmet.jpg")).to eq('')

      # cleanup
      FileUtils.rm_rf('spec/data/helmet')
    end

    # We aren't able to test this because Metis uses Apache's X-Sendfile module
    #   to send the actual bits. So we require a 'production' environment.
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

    it 'will re-download data if the file size does not match' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      # first time it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg/ }

      # At this point there will be a file created locally, but it will be 0 bytes
      #   because no bits are sent. So the file size will not be matched with remote
      #   and calling the command again should result in a re-download.

      # second time it still shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg/ }

      # cleanup
      FileUtils.rm_rf('spec/data/helmet')
    end
  end

  describe MetisShell::Mkdir do
    it 'warns if the path folder exists' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      expect(Metis::Folder.count).to eq(1)

      expect_error("metis://athena/armor", "mkdir", "helmet") { /Folder already exists/ }

      expect(Metis::Folder.count).to eq(1)
    end

    it 'creates a folder hierarchy like mkdir -p' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      expect(Metis::Folder.count).to eq(1)

      expect_output("metis://athena/armor", "mkdir", "shield/shiny/padded") { // }

      expect(Metis::Folder.count).to eq(4)
    end
  end

  describe MetisShell::Rm do
    it 'warns if the path does not exist' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(1)

      expect_error("metis://athena/armor", "rm", "helmet/helmet2.jpg") { /Invalid path/ }

      expect(Metis::File.count).to eq(1)
    end

    it 'warns if the parent folder does not exist' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(1)

      expect_error("metis://athena/armor", "rm", "helmet2/helmet.jpg") { /Invalid folder/ }

      expect(Metis::File.count).to eq(1)
    end

    it 'removes a file' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(1)

      expect_output("metis://athena/armor", "rm", "helmet/helmet.jpg") { // }

      expect(Metis::File.count).to eq(0)
    end

    it 'does nothing if the source folder does not exist' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)

      expect(Metis::Folder.count).to eq(1)

      expect_error("metis://athena/armor", "rm", "helmet2") { /Invalid path/ }

      expect(Metis::Folder.count).to eq(1)
    end

    it 'removes a folder' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)

      expect(Metis::Folder.count).to eq(1)

      expect_output("metis://athena/armor", "rm", "helmet") { // }

      expect(Metis::Folder.count).to eq(0)
    end
  end

  describe MetisShell::Cp do
    it 'copies a file in the same directory' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(1)

      expect_output("metis://athena/armor/helmet", "cp", "helmet.jpg", "helmet2.jpg") { // }

      expect(Metis::File.count).to eq(2)
    end

    it 'copies a file using a relative destination path' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)
      sub_folder = create_folder('athena', 'sub_helmet', bucket: bucket, folder: helmet_folder)

      expect(Metis::File.count).to eq(1)

      expect_output("metis://athena/armor/helmet/sub_helmet", "cp", "../helmet.jpg", ".") { // }

      expect(Metis::File.count).to eq(2)
    end

    it 'copies a file using a relative parent destination path' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      sub_folder = create_folder('athena', 'sub_helmet', bucket: bucket, folder: helmet_folder)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: sub_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/sub_folder/helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(1)

      expect_output("metis://athena/armor/helmet/sub_helmet", "cp", "helmet.jpg", "../") { // }

      expect(Metis::File.count).to eq(2)
    end

    it 'copies a file to a new bucket using metis path' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'editor', owner: 'metis')

      expect(Metis::File.count).to eq(1)

      expect_output("metis://athena/armor/helmet", "cp", "helmet.jpg", "metis://athena/sundry/helmet.jpg") { // }

      expect(Metis::File.count).to eq(2)
    end

    it 'copies a file from a different bucket using metis path' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      sundry_bucket = create( :bucket, project_name: 'athena', name: 'sundry', access: 'editor', owner: 'metis')

      expect(Metis::File.count).to eq(1)

      expect_output("metis://athena/sundry", "cp", "metis://athena/armor/helmet/helmet.jpg", "helmet.jpg") { // }

      expect(Metis::File.count).to eq(2)
    end
  end
end
