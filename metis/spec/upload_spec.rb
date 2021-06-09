require 'digest/md5'

describe UploadController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    default_bucket('athena')

    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}" if set_metis_uid
  end

  let(:set_metis_uid) { true }

  after(:each) do
    stubs.clear
  end


  def upload_path(project_name,file_name)
    "/#{project_name}/upload/files/#{file_name}"
  end

  context '#authorize' do
    it 'should authorize an upload' do
      params = {
        file_path: 'wisdom.txt',
        project_name: 'athena',
        bucket_name: 'files'
      }

      # we use our token to authorize an upload
      token_header(:editor)
      json_post('/authorize/upload', params)

      # we expect an authorization url in return
      uri = URI(json_body[:url])
      hmac_params = Rack::Utils.parse_nested_query(uri.query)

      expect(last_response.status).to eq(200)
      expect(uri.path).to eq("/#{params[:project_name]}/upload/files/#{params[:file_path]}")
      expect(hmac_params['X-Etna-Id']).to eq('metis')

      expect(hmac_params['X-Etna-Email']).to eq('metis@olympus.org')
      expect(hmac_params['X-Etna-Name']).to eq('Metis')
      expect(hmac_params['X-Etna-Headers']).to eq('email,name')
      upload = Metis::Upload.first
      expect(upload).to be_nil  # No Uploads now until upload_start
    end

    it 'should re-use uploads' do
      # there is already an upload in place for our metis_uid
      upload = create_upload( 'athena', 'wisdom.txt', @metis_uid,
        current_byte_position: 10,
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 'abcdef'
      )

      params = {
        file_path: 'wisdom.txt',
        project_name: 'athena',
        bucket_name: 'files'
      }

      # we use our token to authorize an upload
      token_header(:editor)
      json_post('/authorize/upload', params)

      # the upload request is okay despite pre-existing
      expect(last_response.status).to eq(200)

      # we expect an authorization url in return
      uri = URI(json_body[:url])
      expect(uri.path).to eq("/#{params[:project_name]}/upload/files/#{params[:file_path]}")

      hmac_params = Rack::Utils.parse_nested_query(uri.query)
      expect(hmac_params['X-Etna-Id']).to eq('metis')

      # We expect no new upload to be created
      expect(Metis::Upload.count).to eq(1)
    end

    it 'does not authorize poorly-named files' do
      params = {
        project_name: 'athena',
        bucket_name: 'files',
        user_email: 'metis@ucsf.edu'
      }

      # we use our token to authorize an upload
      token_header(:editor)

      # here are some good and bad examples to test
      {
        422 => ['"wisdom".txt''?;wisdom.txt',"wisdom.txt\n\r"],
        200 => [ 'wisdom [1](2){3}.txt' ]
      }.each do |status, examples|
        # we try each example and see if it returns the
        # appropriate status
        examples.each do |example|
          json_post('/authorize/upload', params.merge(file_path: example))
          expect(last_response.status).to eq(status)
        end
      end
    end

    it 'does not authorize read-only files' do
      # there is an existing read-only file
      file = create_file('athena', 'wisdom.txt', WISDOM,
        author: 'metis@tartarus.org',
        read_only: true
      )

      params = {
        file_path: 'wisdom.txt',
        bucket_name: 'files',
        project_name: 'athena'
      }

      # we use our token to authorize an upload
      token_header(:editor)
      json_post('/authorize/upload', params)

      # we expect to be forbidden from uploading
      expect(last_response.status).to eq(403)
      upload = Metis::Upload.first
      expect(upload).to be_nil
    end

    it 'does not authorize files with no parent folder' do
      # there is an existing read-only file
      file = create_file('athena', 'blueprints/wisdom.txt', WISDOM,
        author: 'metis@tartarus.org',
        read_only: true
      )

      params = {
        file_path: 'blueprints/wisdom.txt',
        bucket_name: 'files',
        project_name: 'athena'
      }

      # we use our token to authorize an upload
      token_header(:editor)
      json_post('/authorize/upload', params)

      # we expect to be forbidden from uploading
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid folder: "blueprints"')
      expect(Metis::Upload.count).to eq(0)
    end

    it 'does not authorize files with a read-only parent folder' do
      blueprints_folder = create_folder('athena', 'blueprints', read_only: true)

      params = {
        file_path: 'blueprints/wisdom.txt',
        bucket_name: 'files',
        project_name: 'athena'
      }
      # we use our token to authorize an upload
      token_header(:editor)
      json_post('/authorize/upload', params)

      # we expect to be forbidden from uploading
      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Folder is read-only')
      expect(Metis::Upload.count).to eq(0)
    end

    it 'does not authorize overwriting an existing folder' do
      blueprints_folder = create_folder('athena', 'blueprints')

      params = {
        file_path: 'blueprints',
        bucket_name: 'files',
        project_name: 'athena'
      }
      # we use our token to authorize an upload
      token_header(:editor)
      json_post('/authorize/upload', params)

      # we expect to be forbidden from uploading
      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot overwrite existing folder')
      expect(Metis::Upload.count).to eq(0)
    end
  end

  context '#upload_start' do
    it 'should start an upload' do
      # we expect the appropriate records to have been created
      upload = create_upload( 'athena', 'wisdom.txt', @metis_uid, file_size: WISDOM.length)

      # we post to the upload path with hmac authorization
      hmac_header
      json_post(
        upload_path('athena', 'wisdom.txt'),
        file_size: WISDOM.length,
        action: 'start',
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # the file has not been made yet
      expect(Metis::File.count).to eq(0)

      # we expect the upload as a json object
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(
        current_byte_position: 0,
        file_name: 'wisdom.txt',
        project_name: 'athena',
        author: 'metis@olympus.org|Metis',
        next_blob_hash: '10',
        next_blob_size: 10
      )
    end

    it 'should resume an existing upload' do
      # there is already an upload in place for our metis_uid
      file = create_file('athena', 'wisdom.txt', WISDOM)
      upload = create_upload( 'athena', 'wisdom.txt', @metis_uid,
        current_byte_position: 10,
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 'abcdef'
      )

      # we attempt to initiate a new upload
      hmac_header
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 'defabc'
      )

      # we receive back the status of the existing upload
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 'abcdef',
        author: 'metis@olympus.org|Metis',
        project_name: 'athena',
        file_name: 'wisdom.txt'
      )
    end

    describe 'when a metis_uid is not set' do
      let(:set_metis_uid) { false }

      it 'gracefully handles a lack of a metis_uid' do
        # we attempt to initiate a new upload
        hmac_header
        json_post(
          upload_path('athena', 'wisdom.txt'),
          action: 'start',
          file_size: WISDOM.length,
          next_blob_size: 10,
          next_blob_hash: 'defabc'
        )

        expect(last_response.status).to eq(422)
        expect(json_body).to eq({
          error: "metis_uid not set, did you forget to configure metis_uid to METIS_TEST_UID?"
        })
      end
    end

    context 'when reset=true and the upload exists' do
      it 'should reset the position of the upload' do
        file = create_file('athena', 'wisdom.txt', WISDOM)
        upload = create_upload( 'athena', 'wisdom.txt', @metis_uid,
                                current_byte_position: 10,
                                file_size: WISDOM.length,
                                next_blob_size: 10,
                                next_blob_hash: 'abcdef'
        )

        # we attempt to initiate a new upload
        hmac_header
        json_post(
            upload_path('athena', 'wisdom.txt'),
            action: 'start',
            file_size: WISDOM.length,
            next_blob_size: 3,
            next_blob_hash: 'defabc',
            reset: true,
        )

        # we receive back the status of the new upload
        expect(last_response.status).to eq(200)
        expect(json_body).to eq(
                                 current_byte_position: 0,
                                 next_blob_size: 3,
                                 next_blob_hash: 'defabc',
                                 author: 'metis|metis',
                                 project_name: 'athena',
                                 file_name: 'wisdom.txt'
                             )
      end
    end

    context 'when the upload exists and the filesize differs' do
      it 'should reset the position of the upload' do
        file = create_file('athena', 'wisdom.txt', WISDOM)
        upload = create_upload( 'athena', 'wisdom.txt', @metis_uid,
                                current_byte_position: 10,
                                file_size: WISDOM.length,
                                next_blob_size: 10,
                                next_blob_hash: 'abcdef'
        )

        # we attempt to initiate a new upload
        hmac_header
        json_post(
            upload_path('athena', 'wisdom.txt'),
            action: 'start',
            file_size: WISDOM.length + 3,
            next_blob_size: 3,
            next_blob_hash: 'defabc',
          )

        # we receive back the status of the new upload
        expect(last_response.status).to eq(200)
        expect(json_body).to eq(
                                 current_byte_position: 0,
                                 next_blob_size: 3,
                                 next_blob_hash: 'defabc',
                                 author: 'metis|metis',
                                 project_name: 'athena',
                                 file_name: 'wisdom.txt'
                             )
      end
    end

    it 'should create a new upload on start, if does not exist' do
      file = create_file('athena', 'wisdom.txt', WISDOM)

      # we create an upload, but the metis_uid is different from ours
      upload = create_upload(
        'athena', 'wisdom.txt', @metis_uid.reverse,
        file_size: WISDOM.length,
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      expect(Metis::Upload.count).to eq(1)

      # we attempt to post to the path
      hmac_header
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we get a different upload back
      expect(last_response.status).to eq(200)
      expect(Metis::Upload.count).to eq(2)
    end

    it 'should create a new upload on start using hmac email and name, if does not exist' do
      file = create_file('athena', 'wisdom.txt', WISDOM)

      # we create an upload, but the metis_uid is different from ours
      upload = create_upload(
        'athena', 'wisdom.txt', @metis_uid.reverse,
        file_size: WISDOM.length,
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      expect(Metis::Upload.count).to eq(1)

      # we attempt to post to the path
      hmac_header(params={
        email: 'athena@olympus.org',
        name: 'Athena Pallas'
      })
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we get a different upload back
      expect(last_response.status).to eq(200)
      expect(Metis::Upload.count).to eq(2)
      expect(Metis::Upload.last.author).to eq('athena@olympus.org|Athena Pallas')
    end

    it 'should reject a request without a valid HMAC' do
      file = create_file('athena', 'wisdom.txt', WISDOM)

      # we create an upload, but the metis_uid is different from ours
      upload = create_upload(
        'athena', 'wisdom.txt', @metis_uid.reverse,
        file_size: WISDOM.length,
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      expect(Metis::Upload.count).to eq(1)

      # we attempt to post to the path
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we get a different upload back
      expect(last_response.status).to eq(401)
      expect(Metis::Upload.count).to eq(1)
    end
  end

  context '#upload_blob' do
    it 'should add data to an existing blob' do
      # the existing data
      partial = WISDOM[0..9]
      # the next piece of data
      next_blob = WISDOM[10..19]

      # there is a file record in place
      file = create_file('athena', 'wisdom.txt', WISDOM)

      # the upload expects the contents of next_blob
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,

        file_size: WISDOM.length,
        current_byte_position: partial.length,
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # the current uploaded file-on-disk
      partial_file = stubs.create_partial(upload, partial, @metis_uid)

      # the file-on-disk the client will post
      wisdom_blob_file = stubs.create_data('athena', 'wisdom_blob', next_blob)

      # post the new blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 10,
        next_blob_hash: 10,
        current_byte_position: partial.length,
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      expect(last_response.status).to eq(200)

      expect(File.read(partial_file)).to eq(WISDOM[0..19])
    end

    it 'does not add data with an incorrect hash' do
      # as before, existing and next data
      partial = WISDOM[0..9]
      next_blob = WISDOM[10..19]

      # we create the blob with the wrong contents
      wisdom_blob_file = stubs.create_data('athena', 'wisdom_blob', next_blob.reverse)

      file = create_file('athena', 'wisdom.txt', WISDOM)

      # the upload expects the correct contents
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: partial.length,
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      partial_file = stubs.create_partial(upload, partial, @metis_uid)

      # we post the wrong blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 10,
        next_blob_hash: 10,
        current_byte_position: partial.length,
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      # the server responds with a client error
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Blob integrity failed')
    end
  end

  context '#upload_blob completion' do
    before(:each) do
      # the next blob completes the data
      @next_blob = WISDOM[20..-1]
      @wisdom_blob_file = stubs.create_data('athena', 'wisdom_blob', @next_blob)

      @partial = WISDOM[0..19]
    end

    def prep_upload(file_name)
      # the upload expects the final blob
      upload = create_upload('athena', file_name, @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: @partial.length,
        next_blob_size: @next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(@next_blob)
      )

      @partial_file = stubs.create_partial(upload, @partial, @metis_uid)

      return upload
    end

    def complete_upload(file_name)
      hmac_header
      post(
        upload_path('athena', file_name),
        action: 'blob',
        next_blob_size: 0,
        next_blob_hash: '',
        current_byte_position: @partial.length,
        blob_data: Rack::Test::UploadedFile.new(
          @wisdom_blob_file,
          'application/octet-stream'
        )
      )
    end

    def prep_existing_file(file_name, contents, params={})
      @creation_time = Date.today - 10
      Timecop.freeze(@creation_time)
      file = create_file(
        'athena', file_name, contents,
        {
          author: 'metis@tartarus.org'
        }.merge(params)
      )
      Timecop.return
      @current_file = stubs.create_file('athena', 'files', file_name, contents)
    end

    it 'finishes the upload' do
      upload = prep_upload('wisdom.txt')

      # post the blob with no next blob
      complete_upload('wisdom.txt')

      expect(last_response.status).to eq(200)

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      file = Metis::File.first

      # the actual file exists with the correct content
      expect(File.read(file.data_block.location)).to eq(WISDOM)

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # we get the file hash in the upload with a download_url
      expect(json_body[:file][:download_url]).not_to be_nil

      # the author is set
      expect(file.author).to eq('metis@olympus.org|Metis')

      # clean up the file
      File.delete(file.data_block.location)
      Timecop.return
    end

    it 'completes over an existing file' do
      # the next blob completes the data
      upload = prep_upload('wisdom.txt')

      prep_existing_file('wisdom.txt', WISDOM*2)

      # post the blob with no next blob
      complete_upload('wisdom.txt')

      expect(last_response.status).to eq(200)

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      file = Metis::File.first

      # the actual file exists with the correct content
      expect(File.read(file.data_block.location)).to eq(WISDOM)

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # the created_at timestamp is the original date
      expect(file.created_at).to be_within(1).of(@creation_time.to_time)

      # the updated_at timestamp is new
      expect(file.updated_at).to be_within(1).of(Time.now)

      # so is the author
      expect(file.author).to eq('metis@olympus.org|Metis')

      # clean up the file
      File.delete(file.data_block.location)
    end

    it 'refuses to complete over an existing read-only file' do
      # the next blob completes the data
      upload = prep_upload('wisdom.txt')

      prep_existing_file('wisdom.txt', WISDOM*2, read_only: true)

      # post the blob with no next blob
      complete_upload('wisdom.txt')

      expect(last_response.status).to eq(403)

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      file = Metis::File.first

      # the actual file exists with the original content
      expect(File.read(file.data_block.location)).to eq(WISDOM*2)

      # the file metadata is the same
      expect(file.created_at).to be_within(1).of(@creation_time.to_time)
      expect(file.updated_at).to be_within(1).of(@creation_time.to_time)
      expect(file.author).to eq('metis@tartarus.org')

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # clean up the file
      File.delete(file.data_block.location)
    end

    it 'refuses to complete over an existing folder' do
      # the next blob completes the data
      upload = prep_upload('wisdom.txt')

      wisdom_folder = create_folder('athena', 'wisdom.txt')

      # post the blob with no next blob
      complete_upload('wisdom.txt')

      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('Cannot overwrite existing folder')

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      # the folder still exists
      expect(Metis::Folder.count).to eq(1)

      # the file does not
      expect(Metis::File.count).to eq(0)
    end

    it 'sets a folder name when it completes' do
      blueprints_folder = create_folder('athena', 'blueprints')
      stubs.create_folder('athena', 'files', 'blueprints')
      # the next blob completes the data
      upload = prep_upload('blueprints/wisdom.txt')

      # post the blob with no next blob
      complete_upload('blueprints/wisdom.txt')

      expect(last_response.status).to eq(200)

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      # there is a new file
      file = Metis::File.last

      expect(file).not_to be_nil

      # the actual file exists with the original content
      expect(File.read(file.data_block.location)).to eq(WISDOM)

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # the file metadata is updated
      expect(file.created_at).to be_within(1).of(Time.now)
      expect(file.updated_at).to be_within(1).of(Time.now)

      # so is the author
      expect(file.author).to eq('metis@olympus.org|Metis')

      # crucially, the folder is set
      expect(file.folder).to eq(blueprints_folder)

      # clean up the file
      File.delete(file.data_block.location)
    end

    it 'does not write into a read-only folder' do
      blueprints_folder = create_folder('athena', 'blueprints', read_only: true)
      # the next blob completes the data
      prep_upload('blueprints/wisdom.txt')

      # post the blob with no next blob
      complete_upload('blueprints/wisdom.txt')

      expect(last_response.status).to eq(403)

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      expect(Metis::File.count).to eq(0)
    end

    it 'does not write into a missing folder' do
      # the next blob completes the data
      prep_upload('blueprints/wisdom.txt')

      # post the blob with no next blob
      complete_upload('blueprints/wisdom.txt')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid folder: "blueprints"')

      # the partial is destroyed
      expect(File.exists?(@partial_file)).to be_falsy

      expect(Metis::File.count).to eq(0)
    end
  end

  context '#upload_cancel' do
    it 'should add cancel an upload' do
      # there is partial data
      partial = WISDOM[0..9]

      # there is an upload waiting
      upload = create_upload( 'athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: partial.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      partial_file = stubs.create_partial(upload, partial, @metis_uid)

      # we post a cancel request with our hmac url
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'cancel',
      )

      # the partial is deleted
      expect(File.exists?(partial_file)).to be_falsy
      expect(Metis::Upload.count).to eq(0)

      expect(last_response.status).to eq(200)
    end

    it 'deletes the file record if there is no existing file data' do
      # there is partial data
      partial = WISDOM[0..9]

      # there is an upload waiting
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: partial.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      partial_file = stubs.create_partial(upload, partial, @metis_uid)

      # we post a cancel request with our hmac url
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'cancel',
      )

      # the partial is deleted
      expect(File.exists?(partial_file)).to be_falsy
      expect(Metis::Upload.count).to eq(0)
      expect(Metis::File.count).to eq(0)

      expect(last_response.status).to eq(200)
    end

    it 'keeps the file record if there is existing file data' do
      # there is partial data
      partial = WISDOM[0..9]

      # our file has existing data
      file = create_file('athena', 'wisdom.txt', WISDOM)
      current_file = stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      # there is an upload waiting
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: partial.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      partial_file = stubs.create_partial(upload, partial, @metis_uid)

      # we post a cancel request with our hmac url
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'cancel',
      )

      # the partial is deleted
      expect(File.exists?(partial_file)).to be_falsy
      expect(Metis::Upload.count).to eq(0)
      expect(Metis::File.count).to eq(1)

      expect(last_response.status).to eq(200)

      # clean up
      File.delete(file.data_block.location)
    end
  end
end
