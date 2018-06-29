describe DownloadController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#download' do
    before(:each) do
      @tips = "1. Burn the hydra's neck after cutting.\n2. Use a river to clean the stables."

      @location = stub_file('readme_hercules.txt', @tips, :labors)

      default_bucket('labors')
    end

    after(:each) do
      clear_stubs
    end

    it 'downloads a file' do
      create_file('labors', 'readme_hercules.txt', @tips)

      header(*Etna::TestAuth.hmac_header({}))
      get('/labors/download/files/readme_hercules.txt', { })

      expect(last_response.status).to eq(200)

      # normally our web server should catch this header and replace the
      # contents; we can't do that with Rack::Test
      expect(last_response.headers['X-Sendfile']).to eq(@location)
    end

    it 'fails if there is no file data' do
      clear_stubs
      # we make the file record, but we don't stub the actual file
      create_file('labors', 'readme_hercules.txt', @tips)

      header(*Etna::TestAuth.hmac_header({}))
      get('/labors/download/files/readme_hercules.txt', { })

      expect(last_response.status).to eq(404)
    end

    it 'fails if the file does not exist' do
      # we create no file record

      header(*Etna::TestAuth.hmac_header({}))
      get('/labors/download/files/readme_hercules.txt', { })

      expect(last_response.status).to eq(404)
    end

    it 'refuses without hmac auth' do
      get('/labors/download/files/readme_hercules.txt', { })

      expect(last_response.status).to eq(401)
    end
  end
end
