def stub_file file_name, contents
  file_loc = File.expand_path(File.join(Metis.instance.project_dir(:labors),file_name))
  File.open(file_loc, 'w') do |file|
    file.puts contents
  end
  return file_loc
end

def remove_stubs
  file_loc = Metis.instance.project_dir(:labors)

  Dir.glob(File.join(file_loc, '*')).each do |file|
    File.delete(file)
  end
end


describe DownloadController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  after(:each) do
    remove_stubs
  end

  context '#download' do
    it 'downloads a file' do
      tips = "1. Burn the hydra's neck after cutting.\n2. Use a river to clean the stables."

      location = stub_file('readme_hercules.txt', tips)
      create(:file, project_name: 'labors', file_name: 'readme_hercules.txt', original_name: 'readme_hercules.txt', uploader: 'outis', size: tips.length)

      header(*Etna::TestAuth.hmac_header({}))
      get('/labors/download/readme_hercules.txt', { })

      expect(last_response.status).to eq(200)

      # normally our web server should catch this header and replace the
      # contents; we can't do that with Rack::Test
      expect(last_response.headers['X-Sendfile']).to eq(location)
    end

    it 'fails if there is no file data' do
      create(:file, project_name: 'labors', file_name: 'readme_hercules.txt', original_name: 'readme_hercules.txt', uploader: 'outis', size: 0)

      header(*Etna::TestAuth.hmac_header({}))
      get('/labors/download/readme_hercules.txt', { })

      expect(last_response.status).to eq(404)
    end

    it 'fails if the file does not exist' do
      tips = "1. Burn the hydra's neck after cutting.\n2. Use a river to clean the stables."
      location = stub_file('readme_hercules.txt', tips)

      # we create no file record

      header(*Etna::TestAuth.hmac_header({}))
      get('/labors/download/readme_hercules.txt', { })

      expect(last_response.status).to eq(404)
    end

    it 'refuses without hmac auth' do
      tips = "1. Burn the hydra's neck after cutting.\n2. Use a river to clean the stables."

      stub_file('readme_hercules.txt', tips)
      get('/labors/download/readme_hercules.txt', { })

      expect(last_response.status).to eq(401)
    end
  end
end
