require_relative '../lib/commands'

describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#ledger backfilled" do

    it "returns details for a single project " do
      disable_all_ledger_events
      
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
      expect(wisdom_detail[:description]).to eq("Originally for wisdom.txt")
      expect(wisdom_detail[:files]).to eq([])
      
      helmet_detail = response[:vacuum][:details].find { |d| d[:data_block_id] == helmet_datablock.id }
      expect(helmet_detail).not_to be_nil
      expect(helmet_detail[:md5_hash]).to eq(helmet_datablock.md5_hash)
      expect(helmet_detail[:size]).to eq(HELMET.bytesize)
      expect(helmet_detail[:description]).to eq("Originally for helmet.jpg")
      expect(helmet_detail[:files]).to eq([])
    end

    it "returns details for multiple projects" do
      disable_all_ledger_events
      
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
      expect(wisdom_detail[:description]).to eq("Originally for wisdom.txt")
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
      default_bucket('athena')
      default_bucket('labors')
      default_bucket('backup')
      
      @metis_uid = Metis.instance.sign.uid
      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
      
      @user = Etna::User.new(AUTH_USERS[:editor])
    end

    it "returns vacuum stats with include_projects parameter" do
      enable_all_ledger_events
      
      # Create files in athena
      athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
      shared_datablock = athena_file.data_block

      another_athena_file = upload_file_via_api('athena', 'athena2.txt', "some other content")
      another_datablock = another_athena_file.data_block
      
      # Create file in labors with same content (reuses datablock)
      labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
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

      # Should have event counts (for athena + included projects: labors, backup)
      expect(response[:event_counts][:create_datablock]).to eq(5) # athena.txt, athena2.txt, labors.txt, labors2.txt, backup.txt
      expect(response[:event_counts][:link_file_to_datablock]).to eq(5) # all 5 files linked
      expect(response[:event_counts][:resolve_datablock]).to eq(3) # athena.txt, athena2.txt, labors2.txt resolved (labors.txt and backup.txt reused)
      expect(response[:event_counts][:reuse_datablock]).to eq(2) # labors.txt and backup.txt reused athena's datablock
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(1) # athena.txt deleted

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
      
    end
  end

  context "#ledger mixed mode (backfilled + tracked)" do
    before(:each) do
      default_bucket('athena')
      default_bucket('labors')
      
      @metis_uid = Metis.instance.sign.uid
      set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
      
      @user = Etna::User.new(AUTH_USERS[:editor])
    end

    it "correctly filters backfilled and tracked events when both exist" do
      disable_all_ledger_events
      
      # Phase 1: Create and delete files with ledger disabled (will be backfilled)
      untracked_wisdom = upload_file_via_api('athena', 'untracked_wisdom.txt', WISDOM)
      untracked_helmet = upload_file_via_api('athena', 'untracked_helmet.jpg', HELMET)
      
      backfilled_wisdom_datablock_id = untracked_wisdom.data_block_id
      backfilled_helmet_datablock_id = untracked_helmet.data_block_id
      
      token_header(:editor)
      delete("/athena/file/remove/files/untracked_wisdom.txt")
      delete("/athena/file/remove/files/untracked_helmet.jpg")
      expect(last_response.status).to eq(200)
      
      # Run backfill to create SYSTEM_BACKFILL events for orphaned datablocks
      backfill_ledger = Metis::BackfillDataBlockLedger.new
      backfill_ledger.execute(project_name: 'athena', links: true)
      backfill_ledger.execute(orphaned: true)
      
      # Phase 2: Enable ledger and create/delete files (tracked mode)
      # Create one tracked file with the same content as the untracked wisdom datablock
      enable_all_ledger_events
      
      tracked_file_1 = upload_file_via_api('athena', 'tracked_file_1.txt', WISDOM)
      tracked_file_2 = upload_file_via_api('athena', 'tracked_file_2.txt', "Second tracked content")
      
      tracked_datablock_1_id = tracked_file_1.data_block_id
      tracked_datablock_2_id = tracked_file_2.data_block_id
      
      # Only delete the first file
      token_header(:editor)
      delete("/athena/file/remove/files/tracked_file_1.txt")
      expect(last_response.status).to eq(200)
      
      # Query backfilled stats - should ONLY show backfilled events
      token_header(:supereditor)
      get('/api/stats/ledger', backfilled: true)
      
      expect(last_response.status).to eq(200)
      backfilled_response = json_body
      
      # Backfilled stats should only count SYSTEM_BACKFILL events (no tracked events)
      expect(backfilled_response[:event_counts][:create_datablock]).to eq(0)
      expect(backfilled_response[:event_counts][:link_file_to_datablock]).to eq(0)
      expect(backfilled_response[:event_counts][:resolve_datablock]).to eq(0)
      expect(backfilled_response[:event_counts][:unlink_file_from_datablock]).to eq(2) # backfill unlinks for wisdom and helmet
      expect(backfilled_response[:event_counts][:reuse_datablock]).to eq(0)
      expect(backfilled_response[:event_counts][:remove_datablock]).to eq(0)
      
      # Vacuum stats should show 2 backfilled orphaned datablocks
      expect(backfilled_response[:vacuum][:datablocks_can_vacuum]).to eq(2)
      expect(backfilled_response[:vacuum][:space_can_clear]).to eq(WISDOM.bytesize + HELMET.bytesize)
      
      # Verify vacuum details for backfilled datablocks
      expect(backfilled_response[:vacuum][:details].length).to eq(2)
      # Verify wisdom datablock has the tracked file associated with it (since tracked_file_1 reused it)
      wisdom_detail = backfilled_response[:vacuum][:details].find { |d| d[:data_block_id] == backfilled_wisdom_datablock_id }
      expect(wisdom_detail).not_to be_nil
      expect(wisdom_detail[:files]).to eq([{file_path: "tracked_file_1.txt", bucket_name: "files"}])
      
      # Verify helmet datablock has an empty files array (no files are using it)
      helmet_detail = backfilled_response[:vacuum][:details].find { |d| d[:data_block_id] == backfilled_helmet_datablock_id }
      expect(helmet_detail).not_to be_nil
      expect(helmet_detail[:files]).to eq([])
      
      # Query tracked stats for athena - should ONLY show tracked events for athena
      get('/api/stats/ledger', project_name: 'athena')

      expect(last_response.status).to eq(200)
      tracked_response = json_body
      
      # Tracked stats should only count tracked events for athena (no backfilled events)
      # We created 2 tracked files: tracked_file_1 reused wisdom datablock, tracked_file_2 created new
      # Only tracked_file_1 was deleted
      expect(tracked_response[:event_counts][:create_datablock]).to eq(2) # Even though untracked_wisdom and tracked_1 share the same datablock, it is okay for tracked_1 to create a new datablock since it gets resolved later
      expect(tracked_response[:event_counts][:link_file_to_datablock]).to eq(2)
      expect(tracked_response[:event_counts][:reuse_datablock]).to eq(1) # tracked_file_1 reused wisdom datablock
      expect(tracked_response[:event_counts][:resolve_datablock]).to eq(1) # checksums for tracked_file_1 and tracked_file_2
      expect(tracked_response[:event_counts][:unlink_file_from_datablock]).to eq(1) # tracked_file_1 deleted
      expect(tracked_response[:event_counts][:remove_datablock]).to eq(0)

      # Verify tracked vacuum details
      expect(tracked_response[:vacuum][:details].length).to eq(1)
      
      # Verify wisdom datablock shows tracked_file_1.txt as orphaned
      tracked_wisdom_detail = tracked_response[:vacuum][:details].find { |d| d[:data_block_id] == backfilled_wisdom_datablock_id }
      expect(tracked_wisdom_detail).not_to be_nil
      expect(tracked_wisdom_detail[:files]).to eq([{file_path: "tracked_file_1.txt", bucket_name: "files"}])
      
    end

  end
end
