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

    it "returns verbose vacuum stats with file details for legacy datablocks" do
      token_header(:supereditor)
      get('/api/stats/ledger', legacy: true, verbose: true)

      expect(last_response.status).to eq(200)
      response = json_body

      # Should have at least 2 orphaned datablocks (wisdom and helmet)
      expect(response[:datablocks_can_vacuum]).to be >= 2
      expect(response[:details].length).to eq(response[:datablocks_can_vacuum])

      # Check structure of details
      response[:details].each do |detail|
        expect(detail).to have_key(:data_block_id)
        expect(detail).to have_key(:md5_hash)
        expect(detail).to have_key(:size)
        expect(detail).to have_key(:project_name)
        expect(detail).to have_key(:files)
        expect(detail[:files]).to be_a(Array)
        
        # Each file entry should have file_path and bucket_name
        detail[:files].each do |file|
          expect(file).to have_key(:file_path)
          expect(file).to have_key(:bucket_name)
        end
      end
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
      
      # Should have at least 1 orphaned datablock (wisdom)
      expect(response[:vacuum][:datablocks_can_vacuum]).to be >= 1
      expect(response[:vacuum][:space_can_clear]).to be >= WISDOM.bytesize
      
      # Should have event counts (keys are symbols after JSON parsing)
      expect(response[:event_counts]).to be_a(Hash)
      expect(response[:event_counts][:create_datablock]).to eq(1)
      expect(response[:event_counts][:link_file_to_datablock]).to eq(1)
      expect(response[:event_counts][:resolve_datablock]).to eq(1)
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(1)
      expect(response[:event_counts][:reuse_datablock]).to eq(0)
      expect(response[:event_counts][:remove_datablock]).to eq(0)
    end

    it "returns verbose vacuum stats with file details for a project" do
      # Create file via upload API (triggers log_link automatically)
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      
      # Delete the file via API (this will automatically log unlink event)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      
      # Check verbose vacuum stats for the project
      token_header(:supereditor)
      get('/api/stats/ledger', project_name: 'athena', verbose: true)
      
      expect(last_response.status).to eq(200)
      response = json_body
      
      # Should have at least 1 orphaned datablock (wisdom)
      expect(response[:vacuum][:datablocks_can_vacuum]).to be >= 1
      expect(response[:vacuum][:space_can_clear]).to be >= WISDOM.bytesize
      expect(response[:vacuum][:details]).to be_a(Array)
      expect(response[:vacuum][:details].length).to eq(response[:vacuum][:datablocks_can_vacuum])
      
      # Check structure of details
      response[:vacuum][:details].each do |detail|
        expect(detail).to have_key(:data_block_id)
        expect(detail).to have_key(:md5_hash)
        expect(detail).to have_key(:size)
        expect(detail[:project_name]).to eq('athena')
        expect(detail).to have_key(:files)
        expect(detail[:files]).to be_a(Array)
        
        # Each file entry should have file_path and bucket_name
        detail[:files].each do |file|
          expect(file).to have_key(:file_path)
          expect(file).to have_key(:bucket_name)
          expect(file[:file_path]).to include('wisdom.txt')
        end
      end
      
      # Should still have event counts
      expect(response[:event_counts]).to be_a(Hash)
      expect(response[:event_counts][:create_datablock]).to eq(1)
      expect(response[:event_counts][:link_file_to_datablock]).to eq(1)
      expect(response[:event_counts][:resolve_datablock]).to eq(1)
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(1)
    end
  end
end
