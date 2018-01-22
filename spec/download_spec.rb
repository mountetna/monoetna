def stub_file file_name, contents
  file_loc = Metis.instance.project_dir(:labors)
  File.open(File.join(file_loc, file_name), "w") do |file|
    file.puts contents
  end
end

def remove_stubs
  file_loc = Metis.instance.project_dir(:labors)

  Dir.glob(File.join(file_loc, "*")).each do |file|
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

  context "#download" do
    it "downloads a file" do
      tips = "1. Burn the hydra's neck after cutting.\n2. Use a river to clean the stables."
      create(:file, project_name: 'labors', file_name: 'readme_hercules.txt')
      header(*Etna::TestAuth.hmac_header({}))
      stub_file('readme_hercules.txt', tips)
      get('/labors/download/readme_hercules.txt', { })

      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq(tips)
    end

    it "refuses without hmac auth" do
      tips = "1. Burn the hydra's neck after cutting.\n2. Use a river to clean the stables."

      stub_file('readme_hercules.txt', tips)
      get('/labors/download/readme_hercules.txt', { })

      expect(last_response.status).to eq(401)
    end
  end
end
