describe SFTPFileDiscoveryJob do
  include Rack::Test::Methods

  let(:secrets) {
    {
      sftp_host: "localhost",
    }
  }

  let(:config) {
    {
      "sftp_path" => "/path/to/files",
    }
  }

  before do
    configure_etna_yml
    ENV['TOKEN'] = TEST_TOKEN
  end

  it "fetch the last scan timestamp from the database" do
    require 'pry'; binding.pry
    job = SFTPFileDiscoveryJob.new(config, secrets)
    job.fetch_last_scan
  end

end
