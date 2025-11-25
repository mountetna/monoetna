require_relative '../lib/commands'

describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#ledger backfilled" do
    before(:each) do
      # Disable ledger for backfill tests
      ENV['METIS_LEDGER_TRACKED_MODE_ENABLED'] = 'false'
    end

    it "returns return backfilled detauls for a single project " do
      lifecycle_result = athena_backfilled_lifecycle
      wisdom_datablock = Metis::DataBlock[lifecycle_result[:wisdom_data_block_id]]
      helmet_datablock = Metis::DataBlock[lifecycle_result[:helmet_data_block_id]]
      
      token_header(:supereditor)
      get('/api/stats/ledger', backfilled: true)

      expect(last_response.status).to eq(200)
      response = json_body

      # Should have event_counts
      expect(response[:event_counts][:create_datablock]).to eq(0)
      expect(response[:event_counts][:link_file_to_datablock]).to eq(0)
      expect(response[:event_counts][:resolve_datablock]).to eq(0)
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(2)
      expect(response[:event_counts][:reuse_datablock]).to eq(0)
      expect(response[:event_counts][:remove_datablock]).to eq(0)

      # Should have vacuum stats
      expect(response[:vacuum][:datablocks_can_vacuum]).to be >= 2
      expect(response[:vacuum][:space_can_clear]).to be >= (WISDOM.bytesize + HELMET.bytesize)
      expect(response[:vacuum][:details].length).to eq(2)
      
      # Verify the two datablock details
      wisdom_detail = response[:vacuum][:details].find { |d| d[:data_block_id] == wisdom_datablock.id }
      expect(wisdom_detail).not_to be_nil
      expect(wisdom_detail[:md5_hash]).to eq(wisdom_datablock.md5_hash)
      expect(wisdom_detail[:size]).to eq(WISDOM.bytesize)
      expect(wisdom_detail[:description]).to eq("wisdom.txt")
      expect(wisdom_detail[:files]).to eq([])
      
      helmet_detail = response[:vacuum][:details].find { |d| d[:data_block_id] == helmet_datablock.id }
      expect(helmet_detail).not_to be_nil
      expect(helmet_detail[:md5_hash]).to eq(helmet_datablock.md5_hash)
      expect(helmet_detail[:size]).to eq(HELMET.bytesize)
      expect(helmet_detail[:description]).to eq("helmet.jpg")
      expect(helmet_detail[:files]).to eq([])
    end

    it "returns return backfilled detauls for a multi project" do
      datablock_ids = multi_project_backfilled_lifecycle
      wisdom_datablock = Metis::DataBlock[datablock_ids[:wisdom_data_block_id]]
      
      token_header(:supereditor)
      get('/api/stats/ledger', backfilled: true)

      expect(last_response.status).to eq(200)
      response = json_body

      # Should have event_counts
      expect(response[:event_counts][:create_datablock]).to eq(0)
      expect(response[:event_counts][:link_file_to_datablock]).to eq(2)
      expect(response[:event_counts][:resolve_datablock]).to eq(0)
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(1)
      expect(response[:event_counts][:reuse_datablock]).to eq(0)
      expect(response[:event_counts][:remove_datablock]).to eq(0)

      # Should have vacuum stats
      expect(response[:vacuum][:datablocks_can_vacuum]).to eq(1)
      expect(response[:vacuum][:space_can_clear]).to eq(WISDOM.bytesize)
      expect(response[:vacuum][:details].length).to eq(1)
      
      # Verify the wisdom datablock detail (shared across projects, now orphaned)
      wisdom_detail = response[:vacuum][:details].find { |d| d[:data_block_id] == wisdom_datablock.id }
      expect(wisdom_detail).not_to be_nil
      expect(wisdom_detail[:md5_hash]).to eq(wisdom_datablock.md5_hash)
      expect(wisdom_detail[:size]).to eq(WISDOM.bytesize)
      expect(wisdom_detail[:description]).to eq("wisdom.txt")
      expect(wisdom_detail[:files]).to eq([]) # no files are associated with the datablock
    end


    it "only allows supereditors" do
      token_header(:editor)
      get('/api/stats/ledger', backfilled: true)

      expect(last_response.status).to eq(403)
    end
  end

  context "#ledger tracked" do
    before(:each) do
      # Enable ledger for tracked mode tests
      ENV['METIS_LEDGER_TRACKED_MODE_ENABLED'] = 'true'
      default_bucket('athena')
      
      @metis_uid = Metis.instance.sign.uid
      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
      
      @user = Etna::User.new(AUTH_USERS[:editor])
    end

    it "returns vacuum stats with include_projects parameter" do
      # Set up project buckets
      default_bucket('labors')
      default_bucket('backup')
      
      # Create files in athena
      athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
      shared_datablock = athena_file.data_block

      another_athena_file = upload_file_via_api('athena', 'athena2.txt', "some other content")
      another_datablock = another_athena_file.data_block
      
      # Create file in labors with same content (reuses datablock)
      labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
      expect(labors_file.data_block_id).to eq(shared_datablock.id)

      another_labors_file = upload_file_via_api('labors', 'labors2.txt',  "some other content 2")

      # Create a file in backup with same content (reuses datablock)
      backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)
      expect(backup_file.data_block_id).to eq(shared_datablock.id)
      
      # Delete athena file (orphaned for athena, but still used by labors)
      token_header(:editor)
      delete("/athena/file/remove/files/athena.txt")
      expect(last_response.status).to eq(200)
      
      # Check vacuum stats without include_projects (should be empty - still used by labors)
      token_header(:supereditor)
      
      # Check vacuum stats with include_projects=['labors', 'backup'] (should include the datablock)
      get('/api/stats/ledger', project_name: 'athena', include_projects: ['labors', 'backup'])
      expect(last_response.status).to eq(200)
      response = json_body
      
      # Should have 1 orphaned datablock when including labors and backup
      expect(response[:vacuum][:datablocks_can_vacuum]).to be >= 1
      expect(response[:vacuum][:space_can_clear]).to be >= WISDOM.bytesize
      expect(response[:vacuum][:include_projects]).to eq(['labors', 'backup'])
      
      # Should have project_breakdown showing datablocks in athena, labors, and backup
      expect(response[:vacuum][:project_breakdown][:athena]).to eq(1)
      expect(response[:vacuum][:project_breakdown][:labors]).to eq(1)
      expect(response[:vacuum][:project_breakdown][:backup]).to eq(1)
      
      # Check that details don't have project_name
      response[:vacuum][:details].each do |detail|
        expect(detail).not_to have_key(:project_name)
      end
      
      # Should still have event counts
      expect(response[:event_counts][:create_datablock]).to eq(1)
      expect(response[:event_counts][:link_file_to_datablock]).to eq(1)
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(1)
    end
  end
end
