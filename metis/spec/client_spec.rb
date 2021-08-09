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

    token = Base64.strict_encode64(
      { email: 'metis@olympus.org', name: "Metis", perm: 'a:athena', exp: 253371439590 }.to_json
    )
    @long_lived_token = "something.#{token}"
    ENV['TOKEN'] = @long_lived_token
    stub_request(:any, %r!^https://metis.test/!).to_rack(app)

    stub_request(:get, "https://janus.test/refresh_token")
      .to_return({
        status: 200,
        body: @long_lived_token
      })
  end

  describe MetisShell do
    it 'tells user if token is expired' do
      token = Base64.strict_encode64(
        { email: 'metis@olympus.org', name: "Metis", perm: 'a:athena', exp: 1000 }.to_json
      )
      ENV['TOKEN'] = "something.#{token}"
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      expect_output("metis://athena", "ls") { %r!expired!}
    end

    it 'refreshes the token if close to expiring' do
      frozen_time = 1000
      Timecop.freeze(DateTime.strptime(frozen_time.to_s, "%s"))
      token = Base64.strict_encode64(
        { email: 'metis@olympus.org', name: "Metis", perm: 'a:athena', exp: frozen_time + 1000 }.to_json
      )
      ENV['TOKEN'] = "something.#{token}"

      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      expect_output("metis://athena", "ls") { %r!armor/!}

      expect(WebMock).to have_requested(:get, "https://janus.test/refresh_token")

      expect(ENV["TOKEN"]).to eq(@long_lived_token)

      Timecop.return
    end

    it 'does not refresh the token if not close to expiring' do
      frozen_time = 1000
      Timecop.freeze(DateTime.strptime(frozen_time.to_s, "%s"))
      token = Base64.strict_encode64(
        { email: 'metis@olympus.org', name: "Metis", perm: 'a:athena', exp: frozen_time + 100000 }.to_json
      )
      full_future_token = "something.#{token}"
      ENV['TOKEN'] = full_future_token

      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      expect_output("metis://athena", "ls") { %r!armor/!}

      expect(WebMock).not_to have_requested(:get, "https://janus.test/refresh_token")
      expect(ENV["TOKEN"]).to eq(full_future_token)

      Timecop.return
    end
  end

  describe MetisShell::Set do
    before(:each) do
      create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
    end

    it 'exits upon error, if errexit option is true' do
      stub_metis_list_bucket

      with_temp_stdio do |stdin, stdout|
        stdin.write("set -e\n")
        stdin.write("project athena\n")
        stdin.write("cd armor\n")
        stdin.write("cd jabberwocky\n") # <- this should fail and exit, the folder doesn't exist
        stdin.write("put #{stdin.path} .\n")
        stdin.flush

        shell = MetisShell.new("metis:://athena/armor")
        begin
          replace_stdio(stdin.path, stdout.path) do
            shell.run
          end
        rescue Exception => e
          # Catch the exception / exit that will get thrown by the shell
          #  otherwise, this test exits, too...
        end

        # Shell exited before it even tried processing the put command
        expect(stdout.read.include?("put")).to eq(false)

        # There is NO call to authorize the upload
        expect(WebMock).not_to have_requested(:post, /https:\/\/metis.test\/authorize\/upload/).
          with(body: hash_including({
            "project_name": "athena",
            "bucket_name": "armor",
            "file_path": ::File.basename(stdin.path)
          }))
      end
    end

    it 'does not exit upon error, if errexit option is false' do
      stub_metis_list_bucket

      with_temp_stdio do |stdin, stdout|
        stdin.write("project athena\n")
        stdin.write("cd armor\n")
        stdin.write("cd jabberwocky\n") # <- this folder doesn't exist, but shell ignores that and moves on to the next command
        stdin.write("put #{stdin.path} .\n")
        stdin.flush

        shell = MetisShell.new("metis:://athena/armor")
        replace_stdio(stdin.path, stdout.path) do
          shell.run
        end

        # The PUT command is still executed, since no `set -e` called.
        expect(stdout.read.include?("put")).to eq(true)

        # There IS a call to authorize the upload
        expect(WebMock).to have_requested(:post, /https:\/\/metis.test\/authorize\/upload/).
          with(body: hash_including({
            "project_name": "athena",
            "bucket_name": "armor",
            "file_path": ::File.basename(stdin.path)
          }))
      end
    end
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


    it 'lists files directly' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'armor', 'helmet/helmet.jpg', HELMET)

      expect_output("metis://athena/armor", "ls", "helmet/helmet.jpg") { %r!armor/helmet/helmet.jpg!m }
    end

    it 'lists files directly in long format' do
      Timecop.freeze(DateTime.parse("2020-06-17T04:37"))
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'armor', 'helmet/helmet.jpg', HELMET)

      expect_output("metis://athena/armor", "ls", "-l", "helmet/helmet.jpg") {
        "metis 13 Jun 17 04:37 armor/helmet/helmet.jpg\n"
      }
      Timecop.return
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
          "X-Etna-Id": "metis",
          "X-Etna-Headers": "email,name"
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
    before(:each) do
      # cleanup
      FileUtils.rm_rf('spec/data/helmet')
    end

    after(:each) do
      # cleanup
      FileUtils.rm_rf('spec/data/helmet')
    end

    it 'downloads a folder from a bucket' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      stub_metis_download("/athena/download/armor/helmet/helmet.jpg", HELMET)

      # it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg.*k?B\/s/ }
    end

    it 'does not download data if local file size matches' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      stub_metis_download("/athena/download/armor/helmet/helmet.jpg", HELMET)

      # first time it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg.*k?B\/s/ }

      # second time it shows no progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /(?!k?B\/s).*/ }
    end

    it 'will re-download data if the local file size does not match' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      stub_metis_download("/athena/download/armor/helmet/helmet.jpg", HELMET)

      # first time it shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg.*k?B\/s/ }

      # At this point there will be a file created locally, so we'll tweak the
      #   server contents to be a different size.
      helmet_file = create_file('athena', 'helmet.jpg', SHINY_HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', SHINY_HELMET)
      stub_metis_download("/athena/download/armor/helmet/helmet.jpg", SHINY_HELMET)

      # second time it still shows download progress
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg.*k?B\/s/ }
      expect(::File.read("spec/data/helmet/helmet.jpg")).to eq(SHINY_HELMET)
    end

    it 'will re-start a download if it is incomplete' do
      allow_any_instance_of(MetisShell::Get).to receive(:sleep)
      
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket, folder: helmet_folder)
      stubs.create_file('athena', 'files', 'armor/helmet/helmet.jpg', HELMET)

      # Should get Incomplete message
      expect_error("metis://athena/armor", "get", "helmet", "spec/data") { /Incomplete download/ }

      stub_metis_download("/athena/download/armor/helmet/helmet.jpg", HELMET)

      # second time it will finish the download
      expect_output("metis://athena/armor", "get", "helmet", "spec/data") { /helmet.jpg.*k?B\/s/ }
      expect(::File.read("spec/data/helmet/helmet.jpg")).to eq(HELMET)
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

  describe MetisShell::Peek do
    describe 'throws error if you peek at' do
      it 'buckets' do
        bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
        expect_error("metis://athena", "peek", "0", "1", "armor") { %r!Cannot peek at buckets!}
      end
  
      it 'folders' do
        bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
        helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
        helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket)
        expect_error("metis://athena/armor", "peek", "0", "1", "helmet") { %r!Cannot peek at folders! }
      end

      it 'non-existent files' do
        bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
        helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
        helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket)
        expect_error("metis://athena/armor", "peek", "0", "1", "vaporware.jpg") { %r!File does not exist! }
      end
    end
    
    it 'peeks at a file' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket)
      stubs.create_file('athena', 'armor', 'helmet.jpg', HELMET)

      MetisShell.new("metis://athena", "peek", "4", "5", "armor/helmet.jpg").run

      expect(WebMock).to have_requested(:get, /https:\/\/metis.test\/athena\/download/).
        with(headers: {
          "Range": "bytes=4-8"
        })
    end

    it 'can write peek output to a file' do
      bucket = create( :bucket, project_name: 'athena', name: 'armor', access: 'editor', owner: 'metis')
      helmet_folder = create_folder('athena', 'helmet', bucket: bucket)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET, bucket: bucket)
      stubs.create_file('athena', 'armor', 'helmet.jpg', HELMET)

      tmp = Tempfile.new

      stub_request(:get, /download/)
        .to_return({
          status: 200,
          body: "partial-chunk"
        })

      MetisShell.new("metis://athena", "peek", "4", "5", "armor/helmet.jpg", tmp.path).run

      expect(WebMock).to have_requested(:get, /https:\/\/metis.test\/athena\/download/).
        with(headers: {
          "Range": "bytes=4-8"
        })

      expect(tmp.read).to eq("partial-chunk")

      tmp.close
      tmp.unlink
    end
  end
end
