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
      expect(response[:vacuum][:datablocks_ready]).to be >= 2
      expect(response[:vacuum][:space_ready]).to be >= (WISDOM.bytesize + HELMET.bytesize)
      expect(response[:vacuum][:datablocks_blocked]).to eq(0)
      full_metadata = response[:vacuum][:full_metadata]
      expect(full_metadata.length).to eq(2)
      expect(response[:vacuum][:date_distribution].values.sum).to eq(2)
      size_distribution = response[:vacuum][:size_distribution].transform_keys(&:to_s)
      extension_distribution = response[:vacuum][:extension_distribution].transform_keys(&:to_s)
      expect(size_distribution["0B-1MB"]).to eq(2)
      expect(extension_distribution["txt"]).to eq(1)
      expect(extension_distribution["jpg"]).to eq(1)
      
      # Verify the two datablock metadata records
      wisdom_record = full_metadata.find { |d| d[:data_block_id] == wisdom_datablock.id }
      expect(wisdom_record).not_to be_nil
      expect(wisdom_record[:md5_hash]).to eq(wisdom_datablock.md5_hash)
      expect(wisdom_record[:size]).to eq(WISDOM.bytesize)
      expect(wisdom_record[:description]).to eq("Originally for wisdom.txt")
      expect(wisdom_record[:files]).to eq([])
      
      helmet_record = full_metadata.find { |d| d[:data_block_id] == helmet_datablock.id }
      expect(helmet_record).not_to be_nil
      expect(helmet_record[:md5_hash]).to eq(helmet_datablock.md5_hash)
      expect(helmet_record[:size]).to eq(HELMET.bytesize)
      expect(helmet_record[:description]).to eq("Originally for helmet.jpg")
      expect(helmet_record[:files]).to eq([])
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
      expect(response[:vacuum][:datablocks_ready]).to eq(1)
      expect(response[:vacuum][:space_ready]).to eq(WISDOM.bytesize)
      expect(response[:vacuum][:datablocks_blocked]).to eq(0)
      full_metadata = response[:vacuum][:full_metadata]
      expect(full_metadata.length).to eq(1)
      expect(response[:vacuum][:date_distribution].values.sum).to eq(1)
      size_distribution = response[:vacuum][:size_distribution].transform_keys(&:to_s)
      extension_distribution = response[:vacuum][:extension_distribution].transform_keys(&:to_s)
      expect(size_distribution["0B-1MB"]).to eq(1)
      expect(extension_distribution["txt"]).to eq(1)
      
      # Verify the wisdom datablock metadata (shared across projects, now orphaned)
      wisdom_record = full_metadata.find { |d| d[:data_block_id] == wisdom_datablock.id }
      expect(wisdom_record).not_to be_nil
      expect(wisdom_record[:md5_hash]).to eq(wisdom_datablock.md5_hash)
      expect(wisdom_record[:size]).to eq(WISDOM.bytesize)
      expect(wisdom_record[:description]).to eq("Originally for wisdom.txt")
      expect(wisdom_record[:files]).to eq([]) # no files are associated with the datablock
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

    it "returns vacuum stats with global live-file check" do
      enable_all_ledger_events

      athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
      shared_datablock = athena_file.data_block

      another_athena_file = upload_file_via_api('athena', 'athena2.txt', "some other content")
      another_datablock = another_athena_file.data_block

      # labors and backup reuse the shared datablock
      labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
      backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)
      expect(labors_file.data_block_id).to eq(shared_datablock.id)
      expect(backup_file.data_block_id).to eq(shared_datablock.id)

      # Delete athena's file — shared_datablock is orphaned for athena but still live in labors/backup
      token_header(:editor)
      delete("/athena/file/remove/files/athena.txt")
      expect(last_response.status).to eq(200)

      # Also delete athena2 — another_datablock is safe to vacuum (not shared)
      delete("/athena/file/remove/files/athena2.txt")
      expect(last_response.status).to eq(200)

      token_header(:supereditor)
      get('/api/stats/ledger', project_name: 'athena')
      expect(last_response.status).to eq(200)
      response = json_body

      # Event counts are scoped to athena only
      expect(response[:event_counts][:create_datablock]).to eq(2)       # athena.txt, athena2.txt
      expect(response[:event_counts][:link_file_to_datablock]).to eq(2) # both linked
      expect(response[:event_counts][:unlink_file_from_datablock]).to eq(2) # both deleted

      # another_datablock is safe to vacuum (not shared)
      expect(response[:vacuum][:datablocks_ready]).to eq(1)
      expect(response[:vacuum][:space_ready]).to be >= another_datablock.size

      # shared_datablock is blocked — labors and backup still have live files
      expect(response[:vacuum][:datablocks_blocked]).to eq(1)
      expect(response[:vacuum][:space_blocked]).to be >= shared_datablock.size
      expect(response[:vacuum][:blocked_by_project][:labors]).to eq(1)
      expect(response[:vacuum][:blocked_by_project][:backup]).to eq(1)

      # No include_projects in the response
      expect(response[:vacuum]).not_to have_key(:include_projects)

      response[:vacuum][:full_metadata].each do |record|
        expect(record).not_to have_key(:project_name)
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
      
      enable_all_ledger_events

      # Run backfill to create SYSTEM_BACKFILL events for orphaned datablocks
      backfill_ledger = Metis::BackfillDataBlockLedger.new
      backfill_ledger.execute(project_name: 'athena', links: true)
      backfill_ledger.execute(orphaned: true)
      
      # Phase 2: Enable ledger and create/delete files (tracked mode)
      # Create one tracked file with the same content as the untracked wisdom datablock
      
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
      
      # Vacuum stats should show 2 backfilled orphaned datablocks (both ready — no live files)
      expect(backfilled_response[:vacuum][:datablocks_ready]).to eq(2)
      expect(backfilled_response[:vacuum][:space_ready]).to eq(WISDOM.bytesize + HELMET.bytesize)
      expect(backfilled_response[:vacuum][:datablocks_blocked]).to eq(0)
      full_metadata = backfilled_response[:vacuum][:full_metadata]
      expect(full_metadata.length).to eq(2)
      expect(backfilled_response[:vacuum][:date_distribution].values.sum).to eq(2)
      size_distribution = backfilled_response[:vacuum][:size_distribution].transform_keys(&:to_s)
      extension_distribution = backfilled_response[:vacuum][:extension_distribution].transform_keys(&:to_s)
      expect(size_distribution["0B-1MB"]).to eq(2)
      expect(extension_distribution["txt"]).to eq(1)
      expect(extension_distribution["jpg"]).to eq(1)
      
      # Verify backfilled full metadata records
      # Verify wisdom datablock has the tracked file associated with it (since tracked_file_1 reused it)
      wisdom_record = backfilled_response[:vacuum][:full_metadata].find { |d| d[:data_block_id] == backfilled_wisdom_datablock_id }
      expect(wisdom_record).not_to be_nil
      expect(wisdom_record[:files]).to eq([{file_path: "tracked_file_1.txt", bucket_name: "files"}])
      
      # Verify helmet datablock has an empty files array (no files are using it)
      helmet_record = backfilled_response[:vacuum][:full_metadata].find { |d| d[:data_block_id] == backfilled_helmet_datablock_id }
      expect(helmet_record).not_to be_nil
      expect(helmet_record[:files]).to eq([])
      
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

      # Verify tracked vacuum metadata
      expect(tracked_response[:vacuum][:full_metadata].length).to eq(1)
      
      # Verify wisdom datablock shows tracked_file_1.txt as orphaned
      tracked_wisdom_record = tracked_response[:vacuum][:full_metadata].find { |d| d[:data_block_id] == backfilled_wisdom_datablock_id }
      expect(tracked_wisdom_record).not_to be_nil
      expect(tracked_wisdom_record[:files]).to eq([{file_path: "tracked_file_1.txt", bucket_name: "files"}])
      
    end

  end
end
