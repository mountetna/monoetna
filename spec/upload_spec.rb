require 'digest/md5'

describe UploadController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    clear_stubs
  end


  WISDOM=<<EOT
Although they are
only breath, words
which I command
are immortal
EOT

  def upload_path(project_name,file_name)
    "#{project_name}/upload/#{file_name}"
  end

  def setup_file(project_name, file_name, contents, params={})
    create( :file,
      {
        project_name: project_name, file_name: file_name,
        original_name: file_name, uploader: 'metis', size: contents.length
      }.merge(params)
    )
  end

  context '#authorize' do
    it 'should authorize an upload' do
      params = {
        file_name: 'wisdom.txt',
        project_name: 'athena',
        user_email: 'metis@ucsf.edu'
      }

      # we use our token to authorize an upload
      header(*Etna::TestAuth.token_header(
        email: 'metis@ucsf.edu', perm: 'e:athena'
      ))
      json_post('authorize/upload', params)

      # we expect an authorization url in return
      url = last_response.body
      uri = URI(url)
      hmac_params = Rack::Utils.parse_nested_query(uri.query)

      expect(last_response.status).to eq(200)
      expect(uri.path).to eq("/#{params[:project_name]}/upload/#{params[:file_name]}")
      expect(hmac_params['X-Etna-Id']).to eq('metis')
    end
  end

  context '#upload_start' do
    it 'should start an upload' do
      # we post to the upload path with hmac authorization
      header(*Etna::TestAuth.hmac_header({}))
      json_post(
        upload_path('athena', 'wisdom.txt'),
        file_size: WISDOM.length,
        action: 'start',
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we expect the appropriate records to have been created
      upload = Metis::Upload.first
      file = Metis::File.first
      expect(upload).not_to be_nil
      expect(file).not_to be_nil
      expect(file.uploads).to eq([upload])

      # we expect the upload as a json object
      expect(last_response.status).to eq(200)
      json = json_body(last_response.body)
      expect(json).to eq(
        current_byte_position: 0,
        file_name: 'wisdom.txt',
        next_blob_hash: '10',
        next_blob_size: 10,
        project_name: 'athena'
      )
    end

    it 'should resume an existing upload' do
      # there is already an upload in place for our metis_uid
      file = setup_file('athena', 'wisdom.txt', WISDOM)
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        current_byte_position: 10,
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 'abcdef'
      )

      # we attempt to initiate a new upload
      header(*Etna::TestAuth.hmac_header({}))
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 'defabc'
      )

      # we receive back the status of the existing upload
      expect(last_response.status).to eq(200)
      json = json_body(last_response.body)
      expect(json).to eq(
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 'abcdef',
        project_name: 'athena',
        file_name: 'wisdom.txt'
      )
    end

    it 'should not resume someone elses upload' do
      file = setup_file('athena', 'wisdom.txt', WISDOM, size: 20)

      # we create an upload, but the metis_uid is different from ours
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid.reverse,
        file_size: WISDOM.length,
        current_byte_position: 10,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we attempt to post to the path
      header(*Etna::TestAuth.hmac_header({}))
      json_post(
        upload_path('athena', 'wisdom.txt'),
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we are forbidden
      expect(last_response.status).to eq(403)
      json = json_body(last_response.body)
      expect(json[:error]).to eq('Upload in progress')
    end
  end

  context '#upload_blob' do
    it 'should add data to an existing blob' do
      # the existing data
      partial = WISDOM[0..9]
      # the next piece of data
      next_blob = WISDOM[10..19]

      # there is a file record in place
      file = setup_file('athena', 'wisdom.txt', WISDOM)

      # the current uploaded file-on-disk
      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      # the upload expects the contents of next_blob
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # the file-on-disk the client will post
      wisdom_blob_file = stub_file(:wisdom_blob, next_blob)

      # post the new blob
      header(*Etna::TestAuth.hmac_header({}))
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
      wisdom_blob_file = stub_file(:wisdom_blob, next_blob.reverse)

      file = setup_file('athena', 'wisdom.txt', WISDOM)

      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      # the upload expects the correct contents
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # we post the wrong blob
      header(*Etna::TestAuth.hmac_header({}))
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
      json = json_body(last_response.body)
      expect(json[:error]).to eq('Blob integrity failed')
    end

    it 'finishes the upload' do
      # the next blob completes the data
      next_blob = WISDOM[20..-1]
      partial = WISDOM[0..19]

      wisdom_blob_file = stub_file(:wisdom_blob, next_blob)

      file = setup_file('athena', 'wisdom.txt', WISDOM)

      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      # the upload expects the final blob
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # post the blob with no next blob
      header(*Etna::TestAuth.hmac_header({}))
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

      # the actual file exists with the correct content
      expect(File.read(file.location)).to eq(WISDOM)

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
      file = setup_file('athena', 'wisdom.txt', WISDOM)
      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      # there is an upload waiting
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we post a cancel request with our hmac url
      header(*Etna::TestAuth.hmac_header({}))
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
      file = setup_file('athena', 'wisdom.txt', WISDOM)
      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      # there is an upload waiting
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we post a cancel request with our hmac url
      header(*Etna::TestAuth.hmac_header({}))
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
      file = setup_file('athena', 'wisdom.txt', WISDOM)
      current_file = stub_file('wisdom.txt', WISDOM, :athena)
      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      # there is an upload waiting
      upload = create( :upload,
        file: file,
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # we post a cancel request with our hmac url
      header(*Etna::TestAuth.hmac_header({}))
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
