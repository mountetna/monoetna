require_relative '../lib/commands'

describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#ledger legacy" do
    before(:each) do
      # Disable ledger for backfill tests
      ENV['METIS_LEDGER_ENABLED'] = 'false'
      create_delete_and_backfill_files
    end

    it "returns vacuum stats for legacy datablocks" do
      token_header(:supereditor)
      get('/api/stats/ledger', legacy: true)

      expect(last_response.status).to eq(200)
      response = json_body

      expect(response).to be_a(Hash)
      expect(response[:datablocks_can_vacuum]).to be_a(Integer)
      expect(response[:space_can_clear]).to be_a(Integer)

      # Should have at least 2 orphaned datablocks (wisdom and helmet)
      expect(response[:datablocks_can_vacuum]).to be >= 2
      expect(response[:space_can_clear]).to be >= (WISDOM.bytesize + HELMET.bytesize)
    end

    it "only allows supereditors" do
      token_header(:editor)
      get('/api/stats/ledger', legacy: true)

      expect(last_response.status).to eq(403)
    end
  end

  context "#ledger project" do
    before(:each) do
      # Enable ledger for project-based tests
      ENV['METIS_LEDGER_ENABLED'] = 'true'
      default_bucket('athena')
      
      @metis_uid = Metis.instance.sign.uid
      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
      
      @user = Etna::User.new(AUTH_USERS[:editor])
    end

    it "returns vacuum stats for a project" do
      # Create file via upload API (triggers log_link automatically)
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      
      # Delete the file via API (this will automatically log unlink event)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      
      # Check vacuum stats for the project
      token_header(:supereditor)
      get('/api/stats/ledger', project_name: 'athena')
      
      expect(last_response.status).to eq(200)
      response = json_body
      
      expect(response).to be_a(Hash)
      expect(response[:project_name]).to eq('athena')
      expect(response[:vacuum]).to be_a(Hash)
      expect(response[:vacuum][:datablocks_can_vacuum]).to be_a(Integer)
      expect(response[:vacuum][:space_can_clear]).to be_a(Integer)
      
      # Should have at least 1 orphaned datablock (wisdom)
      expect(response[:vacuum][:datablocks_can_vacuum]).to be >= 1
      expect(response[:vacuum][:space_can_clear]).to be >= WISDOM.bytesize
      
      # Should have event counts
      expect(response[:event_counts]).to be_a(Hash)
      expect(response[:event_counts][Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK]).to eq(1)
      expect(response[:event_counts][Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK]).to eq(1)
    end
  end
end
