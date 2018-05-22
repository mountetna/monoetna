require 'digest/md5'

describe UploadController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    default_bucket('athena')

    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    clear_stubs
  end


  def upload_path(project_name,file_name)
    "#{project_name}/upload/files/#{file_name}"
  end

  context '#authorize' do
    it 'should authorize an upload' do
      params = {
        file_name: 'wisdom.txt',
        project_name: 'athena'
      }

      # we use our token to authorize an upload
      auth_header(:editor)
      json_post('authorize/upload', params)

      # we expect an authorization url in return
      url = last_response.body
      uri = URI(url)
      hmac_params = Rack::Utils.parse_nested_query(uri.query)

      expect(last_response.status).to eq(200)
      expect(uri.path).to eq("/#{params[:project_name]}/upload/files/#{params[:file_name]}")
      expect(hmac_params['X-Etna-Id']).to eq('metis')
      upload = Metis::Upload.first
      expect(upload).not_to be_nil
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
        file_name: 'wisdom.txt',
        project_name: 'athena'
      }

      # we use our token to authorize an upload
      auth_header(:editor)
      json_post('authorize/upload', params)

      # the upload request is okay despite pre-existing
      expect(last_response.status).to eq(200)

      # we expect an authorization url in return
      uri = URI(last_response.body)
      expect(uri.path).to eq("/#{params[:project_name]}/upload/files/#{params[:file_name]}")

      hmac_params = Rack::Utils.parse_nested_query(uri.query)
      expect(hmac_params['X-Etna-Id']).to eq('metis')

      # We expect no new upload to be created
      expect(Metis::Upload.count).to eq(1)
    end

    it 'does not authorize poorly-named files' do
      params = {
        project_name: 'athena',
        user_email: 'metis@ucsf.edu'
      }

      # we use our token to authorize an upload
      auth_header(:editor)

      # here are some good and bad examples to test
      {
        422 => ['"wisdom".txt''?;wisdom.txt',"wisdom.txt\n\r"],
        200 => [ 'wisdom [1](2){3}.txt' ]
      }.each do |status, examples|
        # we try each example and see if it returns the
        # appropriate status
        examples.each do |example|
          json_post('authorize/upload', params.merge(file_name: example))
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
        file_name: 'wisdom.txt',
        project_name: 'athena'
      }

      # we use our token to authorize an upload
      auth_header(:editor)
      json_post('authorize/upload', params)

      # we expect to be forbidden from uploading
      expect(last_response.status).to eq(403)
      upload = Metis::Upload.first
      expect(upload).to be_nil
    end
  end

  context '#upload_start' do
    it 'should start an upload' do
      # we expect the appropriate records to have been created
      upload = create_upload( 'athena', 'wisdom.txt', @metis_uid )

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
        author: 'metis@olympus.org|Metis ',
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
        author: 'metis@olympus.org|Metis ',
        project_name: 'athena',
        file_name: 'wisdom.txt'
      )
    end

    it 'should not resume someone elses upload' do
      file = create_file('athena', 'wisdom.txt', WISDOM)

      # we create an upload, but the metis_uid is different from ours
      upload = create_upload(
        'athena', 'wisdom.txt', @metis_uid.reverse,
        file_size: WISDOM.length,
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we attempt to post to the path
      hmac_header
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we are forbidden
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('No matching upload!')
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

      # the current uploaded file-on-disk
      partial_file = stub_partial(file.file_name, partial, :athena)


      # the upload expects the contents of next_blob
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # the file-on-disk the client will post
      wisdom_blob_file = stub_data('wisdom_blob', next_blob, :athena)

      # post the new blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 10,
        next_blob_hash: 10,
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
      wisdom_blob_file = stub_data('wisdom_blob', next_blob.reverse, :athena)

      file = create_file('athena', 'wisdom.txt', WISDOM)

      partial_file = stub_partial(file.file_name, partial, :athena)

      # the upload expects the correct contents
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # we post the wrong blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 10,
        next_blob_hash: 10,
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      # the server responds with a client error
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Blob integrity failed')
    end

    it 'finishes the upload' do
      # the next blob completes the data
      next_blob = WISDOM[20..-1]
      partial = WISDOM[0..19]

      wisdom_blob_file = stub_data('wisdom_blob', next_blob, :athena)

      partial_file = stub_partial('wisdom.txt', partial, :athena)

      # the upload expects the final blob
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # post the blob with no next blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 0,
        next_blob_hash: '',
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      expect(last_response.status).to eq(200)

      # the partial is destroyed
      expect(File.exists?(partial_file)).to be_falsy

      file = Metis::File.first

      # the actual file exists with the correct content
      expect(File.read(file.location)).to eq(WISDOM)

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # the author is set
      expect(file.author).to eq('metis@olympus.org|Metis ')

      # clean up the file
      File.delete(file.location)
      Timecop.return
    end

    it 'completes over an existing file' do
      # the next blob completes the data
      next_blob = WISDOM[20..-1]
      partial = WISDOM[0..19]

      wisdom_blob_file = stub_data('wisdom_blob', next_blob, :athena)

      partial_file = stub_partial('wisdom.txt', partial, :athena)

      creation_time = Date.today - 10
      Timecop.freeze(creation_time)
      file = create_file('athena', 'wisdom.txt', WISDOM*2,
        author: 'metis@tartarus.org'
      )
      Timecop.return

      current_file = stub_file('wisdom.txt', WISDOM*2, :athena)

      # the upload expects the final blob
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # post the blob with no next blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 0,
        next_blob_hash: '',
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      expect(last_response.status).to eq(200)

      # the partial is destroyed
      expect(File.exists?(partial_file)).to be_falsy

      file = Metis::File.first

      # the actual file exists with the correct content
      expect(File.read(file.location)).to eq(WISDOM)

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # the created_at timestamp is the original date
      expect(file.created_at).to be_within(1).of(creation_time.to_time)

      # the updated_at timestamp is new
      expect(file.updated_at).to be_within(1).of(Time.now)

      # so is the author
      expect(file.author).to eq('metis@olympus.org|Metis ')

      # clean up the file
      File.delete(file.location)
    end

    it 'refuses to complete over an existing read-only file' do
      # the next blob completes the data
      next_blob = WISDOM[20..-1]
      partial = WISDOM[0..19]

      wisdom_blob_file = stub_data('wisdom_blob', next_blob, :athena)

      partial_file = stub_partial('wisdom.txt', partial, :athena)

      creation_time = Date.today - 10
      Timecop.freeze(creation_time)
      file = create_file('athena', 'wisdom.txt', WISDOM*2,
        author: 'metis@tartarus.org',
        read_only: true
      )
      Timecop.return
      current_file = stub_file('wisdom.txt', WISDOM*2, :athena)

      # the upload expects the final blob
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # post the blob with no next blob
      hmac_header
      post(
        upload_path('athena', 'wisdom.txt'),
        action: 'blob',
        next_blob_size: 0,
        next_blob_hash: '',
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      expect(last_response.status).to eq(403)

      # the partial is destroyed
      expect(File.exists?(partial_file)).to be_falsy

      file = Metis::File.first

      # the actual file exists with the original content
      expect(File.read(file.location)).to eq(WISDOM*2)

      # the file metadata is the same
      expect(file.created_at).to be_within(1).of(creation_time.to_time)
      expect(file.updated_at).to be_within(1).of(creation_time.to_time)
      expect(file.author).to eq('metis@tartarus.org')

      # the file thinks it has data
      expect(file.has_data?).to be_truthy

      # clean up the file
      File.delete(file.location)
    end
  end

  context '#upload_cancel' do
    it 'should add cancel an upload' do
      # there is partial data
      partial = WISDOM[0..9]
      partial_file = stub_partial('wisdom.txt', partial, :athena)

      # there is an upload waiting
      upload = create_upload( 'athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: 10,
        next_blob_hash: 10
      )

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

      # our file has no existing data
      partial_file = stub_partial('wisdom.txt', partial, :athena)

      # there is an upload waiting
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: 10,
        next_blob_hash: 10
      )

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
      current_file = stub_file('wisdom.txt', WISDOM, :athena)
      partial_file = stub_partial(file.file_name, partial, :athena)

      # there is an upload waiting
      upload = create_upload('athena', 'wisdom.txt', @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: 10,
        next_blob_hash: 10
      )

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
      File.delete(file.location)
    end
  end
end
