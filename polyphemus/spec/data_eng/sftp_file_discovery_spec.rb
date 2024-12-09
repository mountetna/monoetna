describe SFTPFileDiscoveryJob do
  include Rack::Test::Methods

  let(:secrets) {
    {
     "sftp_host" => "sftp",
      "sftp_user" => "user",
      "sftp_password" => "password",
      "sftp_port" => "22",
    }
  }

  let(:config) {
    {
      "sftp_path" => "/home/user/uploads",
    }
  }

  before do
    ENV['TOKEN'] = TEST_TOKEN
    ENV['ARGO_WORKFLOW_ID'] = "1234567890"
    # Setup sftp client
  end

  context 'first run' do

    before do
      # clear the db
      job = SFTPFileDiscoveryJob.new(config, secrets)
      job.process 
    end

    it "records the last scan timestamp in the db" do
      # check the db
    end

    it 'records the number of files to update in the db' do
    end

    it 'writes the files to update to a csv' do
    end 


  end

  context 'subsequent runs' do

    before do
      # set the db
    end

    it "fetches the last scan timestamp from the database" do
      job = SFTPFileDiscoveryJob.new(config, secrets)
      job.fetch_last_scan
    end
  end

end
