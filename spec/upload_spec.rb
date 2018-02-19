require 'digest/md5'

describe UploadController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  WISDOM=<<EOT
Although they are
only breath, words
which I command
are immortal
EOT

  context '#authorize' do
    it 'should authorize an upload' do
      params = {
        file_name: 'wisdom.txt',
        project_name: 'athena',
        user_email: 'metis@ucsf.edu'
      }

      header(
        *Etna::TestAuth.token_header(
          email: 'metis@ucsf.edu', perm: 'e:athena'
        )
      )

      json_post('authorize/upload', params)

      url = last_response.body
      uri = URI(url)

      hmac_params = Rack::Utils.parse_nested_query(uri.query)

      expect(last_response.status).to eq(200)
      expect(uri.path).to eq('/upload')
      expect(hmac_params['X-Etna-Id']).to eq('metis')
      expect(hmac_params['X-Etna-Project-Name']).to eq(params[:project_name])
      expect(hmac_params['X-Etna-File-Name']).to eq(params[:file_name])
    end
  end

  context '#upload_start' do
    before(:each) do
      @metis_uid = Metis.instance.sign.uid

      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
    end

    it 'should start an upload' do
      header(*Etna::TestAuth.hmac_header({}))
      json_post(
        'upload',
        project_name: 'athena',
        file_name: 'wisdom.txt',
        file_size: WISDOM.length,
        action: 'start',
        next_blob_size: 10,
        next_blob_hash: 10
      )

      expect(last_response.status).to eq(200)
      json = json_body(last_response.body)
      upload = Metis::Upload.first
      file = Metis::File.first

      expect(upload).not_to be_nil
      expect(file).not_to be_nil
      expect(file.uploads).to eq([upload])

      expect(json[:file_name]).to eq('wisdom.txt')
      expect(json[:project_name]).to eq('athena')
      expect(json[:status]).to eq('initialized')
    end

    it 'should resume an existing upload' do
      # before we start, we create the file and the upload
      file = create( :file, 
        project_name: 'athena', file_name: 'wisdom.txt',
        original_name: 'wisdom.txt', uploader: 'metis', size: WISDOM.length
      )
      upload = create( :upload,
        file: file, status: 'queued',
        metis_uid: @metis_uid,
        current_byte_position: 10,
        current_blob_size: 10,
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # now we post to the same place
      header(*Etna::TestAuth.hmac_header(
        project_name: 'athena',
        file_name: 'wisdom.txt'
      ))
      json_post(
        'upload',
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      expect(last_response.status).to eq(200)

      json = json_body(last_response.body)
      expect(json[:file_name]).to eq('wisdom.txt')
      expect(json[:project_name]).to eq('athena')
      expect(json[:status]).to eq('queued')
      expect(json[:current_byte_position]).to eq(10)
    end

    it 'should not resume someone elses upload' do
      # before we start, we create the file and the upload
      file = create( :file, 
        project_name: 'athena', file_name: 'wisdom.txt',
        size: WISDOM.length,
        original_name: 'wisdom.txt', uploader: 'metis', size: 20
      )
      upload = create( :upload,
        file: file, status: 'queued',
        metis_uid: @metis_uid.reverse,
        file_size: WISDOM.length,
        current_byte_position: 10,
        current_blob_size: 10,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      # now we post to the same place
      header(*Etna::TestAuth.hmac_header(
        project_name: 'athena',
        file_name: 'wisdom.txt'
      ))
      json_post(
        'upload',
        action: 'start',
        file_size: WISDOM.length,
        next_blob_size: 10,
        next_blob_hash: 10
      )

      expect(last_response.status).to eq(422)

      json = json_body(last_response.body)
      expect(json[:error]).to eq('Upload in progress')
    end
  end

  context '#upload_blob' do
    before(:each) do
      @metis_uid = Metis.instance.sign.uid

      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
    end

    after(:each) do
      clear_stubs
    end

    it 'should add data to an existing blob' do
      # before we start, we create the file and the upload
      next_blob = WISDOM[10..19]
      partial = WISDOM[0..9]

      wisdom_blob_file = stub_file(:wisdom_blob, next_blob)

      file = create(:file,
        project_name: 'athena', file_name: 'wisdom.txt',
        size: WISDOM.length,
        original_name: 'wisdom.txt', uploader: 'metis'
      )

      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      upload = create( :upload,
        file: file, status: 'queued',
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        current_blob_size: File.size(wisdom_blob_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # now we post to the same place
      header(*Etna::TestAuth.hmac_header(
        project_name: 'athena',
        file_name: 'wisdom.txt'
      ))

      post(
        '/upload',
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
      # before we start, we create the file and the upload
      next_blob = WISDOM[10..19]
      partial = WISDOM[0..9]

      wisdom_blob_file = stub_file(:wisdom_blob, next_blob.reverse)

      file = create( :file, 
        project_name: 'athena', file_name: 'wisdom.txt',
        size: WISDOM.length,
        original_name: 'wisdom.txt', uploader: 'metis'
      )

      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      upload = create( :upload,
        file: file, status: 'queued',
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        current_blob_size: File.size(wisdom_blob_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # now we post to the same place
      header(*Etna::TestAuth.hmac_header(
        project_name: 'athena',
        file_name: 'wisdom.txt'
      ))

      post(
        '/upload',
        action: 'blob',
        next_blob_size: 10,
        next_blob_hash: 10,
        blob_data: Rack::Test::UploadedFile.new(
          wisdom_blob_file,
          'application/octet-stream'
        )
      )

      expect(last_response.status).to eq(422)
      json = json_body(last_response.body)
      expect(json[:error]).to eq('Blob integrity failed')
    end

    it 'finishes the upload' do
      # before we start, we create the file and the upload
      next_blob = WISDOM[20..-1]
      partial = WISDOM[0..19]

      wisdom_blob_file = stub_file(:wisdom_blob, next_blob)

      file = create( :file, 
        project_name: 'athena', file_name: 'wisdom.txt',
        size: WISDOM.length,
        original_name: 'wisdom.txt', uploader: 'metis'
      )

      partial_file = stub_file(
        "uploads/#{@metis_uid}-#{file.file_name}", partial, :athena
      )

      upload = create( :upload,
        file: file, status: 'queued',
        metis_uid: @metis_uid,
        file_size: WISDOM.length,
        current_byte_position: File.size(partial_file),
        current_blob_size: File.size(wisdom_blob_file),
        next_blob_size: next_blob.length,
        next_blob_hash: Digest::MD5.hexdigest(next_blob)
      )

      # now we post to the same place
      header(*Etna::TestAuth.hmac_header(project_name: 'athena', file_name: 'wisdom.txt'))

      post(
        '/upload',
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
    end
  end
end
