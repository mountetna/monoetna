require_relative 'spec_helper'

describe Metis::DataBlockLedger do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    # Our projects
    default_bucket('athena')
    default_bucket('labors')
    default_bucket('backup')
    
    @metis_uid = Metis.instance.sign.uid
    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
    
    stub_event_log(:athena)
    stub_event_log(:labors)
    stub_event_log(:backup)
  end

  after(:each) do
    stubs.clear
  end

  describe '.find_orphaned_datablocks_backfilled' do
    describe 'single project' do
      it 'returns orphaned datablocks' do
        disable_all_ledger_events
        
        # Create files
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        helmet_file = upload_file_via_api('athena', 'helmet.jpg', HELMET)
        
        wisdom_datablock = wisdom_file.data_block
        helmet_datablock = helmet_file.data_block
        
        # Delete files via API
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        delete("/athena/file/remove/files/helmet.jpg")
        expect(last_response.status).to eq(200)
        
        # Run backfill to create SYSTEM_BACKFILL unlink events
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(project_name: 'athena', links: true)
        backfill_ledger.execute(orphaned: true)
        
        # Find orphaned backfilled datablocks
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
        
        expect(orphaned.length).to eq(2)
        orphaned_ids = orphaned.map(&:id)
        expect(orphaned_ids).to include(wisdom_datablock.id, helmet_datablock.id)
      end

      it 'does not return datablocks that have been vacuumed' do
        disable_all_ledger_events
        
        # Create a file
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        # Delete file via API
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Run backfill to create SYSTEM_BACKFILL unlink event
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(orphaned: true)
        
        # Vacuum the datablock via backfilled API (marks it as removed and creates REMOVE_DATABLOCK event)
        token_header(:supereditor)
        json_post('/api/vacuum_datablocks/backfilled', { commit: true })
        expect(last_response.status).to eq(200)
        
        # Should not be orphaned because it's been vacuumed (marked as removed)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
        
        expect(orphaned).to be_empty
      end

      it 'does not return datablocks that have been triggered in the tracked mode' do
        disable_all_ledger_events
        
         # Create a untracked file
        untracked_file = upload_file_via_api('athena', 'untracked1.txt', WISDOM)
        token_header(:editor)
        delete("/athena/file/remove/files/untracked1.txt")

        # Enable ledger tracking (allow real ledger events)
        enable_all_ledger_events

        # Create a tracked file
        tracked_file = upload_file_via_api('athena', 'tracked2.txt', HELMET)
        delete("/athena/file/remove/files/tracked2.txt")

        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(project_name: 'athena', links: true)
        backfill_ledger.execute(orphaned: true)

        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
        expect(orphaned.length).to eq(1)
        expect(orphaned[0].id).to eq(untracked_file.data_block.id)
      end
    end

    describe 'cross project' do
      it 'finds orphaned datablocks that are still in use by another project' do
        disable_all_ledger_events
        
        # Create files with the same data block
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
    
        # Use supereditor to have access to all projects
        token_header(:supereditor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        delete("/labors/file/remove/files/labors.txt")
        expect(last_response.status).to eq(200)

        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
    
        # Run backfill should detect no orphaned datablocks
        backfill_ledger.execute(project_name: 'athena', links: true)
        backfill_ledger.execute(project_name: 'labors', links: true)
        backfill_ledger.execute(orphaned: true)

        # Find orphaned datablocks
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(wisdom_file.data_block_id)
      end

      it 'does not find orphaned datablocks that are still in use by another project' do
        disable_all_ledger_events
        
        # Create files with the same data block
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
        backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)
    
        # Use supereditor to have access to all projects
        token_header(:supereditor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
    
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
    
        # Run backfill should detect no orphaned datablocks
        backfill_ledger.execute(project_name: 'athena', links: true)
        backfill_ledger.execute(project_name: 'labors', links: true)
        backfill_ledger.execute(project_name: 'backup', links: true)
        backfill_ledger.execute(orphaned: true)

        # Find orphaned datablocks
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_backfilled
        expect(orphaned.length).to eq(0)
      end
    end
  end

  describe '.find_orphaned_datablocks' do
    describe 'single project' do
      it 'returns orphaned datablocks that were linked and then unlinked' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create a file via API (triggers log_link automatically)
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        # Verify link event was created
        link_event = Metis::DataBlockLedger.where(
          data_block_id: wisdom_datablock.id,
          event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
        ).first
        expect(link_event).to be_present
        
        # Delete the file via API (triggers log_unlink automatically)
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Verify unlink event was created
        unlink_event = Metis::DataBlockLedger.where(
          data_block_id: wisdom_datablock.id,
          event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
        ).first
        expect(unlink_event).to be_present
        
        # Find orphaned datablocks
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(wisdom_datablock.id)
        expect(orphaned.first.md5_hash).to eq(wisdom_datablock.md5_hash)
      end

      it 'returns orphaned datablocks that were reused and then unlinked' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create first file
        wisdom_file1 = upload_file_via_api('athena', 'wisdom1.txt', WISDOM)
        wisdom_datablock = wisdom_file1.data_block
        
        # Create second file with same content (triggers reuse)
        wisdom_file2 = upload_file_via_api('athena', 'wisdom2.txt', WISDOM)
        expect(wisdom_file2.data_block_id).to eq(wisdom_datablock.id)
        
        # Verify reuse event was created
        reuse_event = Metis::DataBlockLedger.where(
          data_block_id: wisdom_datablock.id,
          event_type: Metis::DataBlockLedger::REUSE_DATABLOCK
        ).first
        expect(reuse_event).to be_present
        
        # Delete both files
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom1.txt")
        delete("/athena/file/remove/files/wisdom2.txt")
        expect(last_response.status).to eq(200)
        
        # Find orphaned datablocks
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(wisdom_datablock.id)
      end

      it 'does not return datablocks still referenced by files' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create a file
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        # Create another file with same content (reuses the datablock)
        wisdom_file2 = upload_file_via_api('athena', 'wisdom2.txt', WISDOM)
        expect(wisdom_file2.data_block_id).to eq(wisdom_datablock.id)
        
        # Delete only one file
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Should not be orphaned because wisdom_file2 still references it
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        
        expect(orphaned).to be_empty
      end

      it 'does not return datablocks that have already been vacuumed' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create and delete a file
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Vacuum the datablock (creates REMOVE_DATABLOCK event)
        Metis::DataBlockLedger.log_vacuum(wisdom_datablock, 'athena', nil)
        
        # Should not be orphaned because it's already been vacuumed
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        
        expect(orphaned).to be_empty
      end

      it 'handles datablocks with multiple link and unlink events correctly' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create a file
        wisdom_file1 = upload_file_via_api('athena', 'wisdom1.txt', WISDOM)
        wisdom_datablock = wisdom_file1.data_block
        
        # Delete it
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom1.txt")
        expect(last_response.status).to eq(200)
        
        # Create another file with same content (reuses datablock)
        wisdom_file2 = upload_file_via_api('athena', 'wisdom2.txt', WISDOM)
        expect(wisdom_file2.data_block_id).to eq(wisdom_datablock.id)
        
        # Delete it too
        delete("/athena/file/remove/files/wisdom2.txt")
        expect(last_response.status).to eq(200)
        
        # Should be orphaned now
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(wisdom_datablock.id)
      end

      it 'returns empty array when no orphaned datablocks exist' do
        disable_all_ledger_events
        
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned).to be_empty
      end

      it 'does not return datablocks that have been triggered in backfilled mode' do
        disable_all_ledger_events
        
         # Create a untracked file
        untracked_file = upload_file_via_api('athena', 'untracked1.txt', WISDOM)
        token_header(:editor)
        delete("/athena/file/remove/files/untracked1.txt")

        # Enable ledger tracking (allow real ledger events)
        enable_all_ledger_events

        # Create a tracked file
        tracked_file = upload_file_via_api('athena', 'tracked2.txt', HELMET)
        delete("/athena/file/remove/files/tracked2.txt")

        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(project_name: 'athena', links: true)
        backfill_ledger.execute(orphaned: true)

        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned.length).to eq(1)
        expect(orphaned[0].id).to eq(tracked_file.data_block.id)
      end
    end

    describe 'cross project' do
      it 'does not return datablocks that are still in use by another project' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create a file in athena project
        athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
        shared_datablock = athena_file.data_block
        
        # Create a file in labors project with same content (reuses the datablock)
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
        expect(labors_file.data_block_id).to eq(shared_datablock.id)

        # Delete the athena file (datablock is orphaned for athena)
        token_header(:editor)
        delete("/athena/file/remove/files/athena.txt")
        expect(last_response.status).to eq(200)
        
        # Should not be orphaned for athena because it's still in use by labors project
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned).to be_empty
      end

      it 'returns datablocks orphaned for project when include_projects list contains all using projects' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create file in athena
        athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
        shared_datablock = athena_file.data_block
        
        # Create file in labors with same content (reuses datablock)
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
        expect(labors_file.data_block_id).to eq(shared_datablock.id)
        
        # Delete athena file (orphaned for athena, but still used by labors)
        token_header(:editor)
        delete("/athena/file/remove/files/athena.txt")
        expect(last_response.status).to eq(200)
        
        # With include=['labors']: should be orphaned (only used by athena + labors)
        orphaned_with_include = Metis::DataBlockLedger.find_orphaned_datablocks('athena', include_projects: ['labors'])
        expect(orphaned_with_include.length).to eq(1)
        expect(orphaned_with_include.first.id).to eq(shared_datablock.id)
      end

      it 'does not return datablocks used by projects outside the include list' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create file in athena
        athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
        shared_datablock = athena_file.data_block
        
        # Create file in labors with same content
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
        expect(labors_file.data_block_id).to eq(shared_datablock.id)
        
        # Create file in backup with same content
        backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)
        expect(backup_file.data_block_id).to eq(shared_datablock.id)
        
        # Delete athena file
        token_header(:editor)
        delete("/athena/file/remove/files/athena.txt")
        expect(last_response.status).to eq(200)
        
        # With include=['labors']: should NOT be orphaned (still used by backup, which is not included)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena', include_projects: ['labors'])
        expect(orphaned).to be_empty
      end

      it 'works with multiple projects in include list' do
        # Allow ledger logging for this test
        enable_all_ledger_events
        
        # Create file in athena
        athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
        shared_datablock = athena_file.data_block
        
        # Create file in labors
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
        expect(labors_file.data_block_id).to eq(shared_datablock.id)
        
        # Create file in backup
        backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)
        expect(backup_file.data_block_id).to eq(shared_datablock.id)
        
        # Delete athena file
        token_header(:editor)
        delete("/athena/file/remove/files/athena.txt")
        expect(last_response.status).to eq(200)
        
        # With include=['labors', 'backup']: should be orphaned (only used by athena + labors + backup)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena', include_projects: ['labors', 'backup'])
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(shared_datablock.id)
      end
    end
  end

  describe 'backfiller and tracked' do
    describe 'mixed mode scenarios' do
      it 'backfiller creates link event, then ledger creates unlink event' do
        # Scenario: File exists when backfill runs, then gets deleted with tracking on
        disable_all_ledger_events
        
        # Upload file without tracking
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        # Run backfill to create link event
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(project_name: 'athena', links: true)
        
        # Enable tracking and delete file (creates real-time unlink event)
        enable_all_ledger_events
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Should correctly identify as orphaned
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(wisdom_datablock.id)
      end

      it 'backfiller creates link event, ledger creates duplicate link event, then unlink event' do
        # Scenario: Backfill runs, then tracking turns on and somehow creates duplicate link, then deletes
        disable_all_ledger_events
        
        # Upload file without tracking
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        # Run backfill to create link event
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(project_name: 'athena', links: true)

        enable_all_ledger_events
        # Upload another file with same content (triggers duplicate link event)
        upload_file_via_api('athena', 'wisdom2.txt', WISDOM)
        
        # Verify we have 1 link event and 1 reuse event
        link_count = Metis::DataBlockLedger.where(
          data_block_id: wisdom_datablock.id,
          event_type: [Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK, Metis::DataBlockLedger::REUSE_DATABLOCK]
        ).count
        expect(link_count).to eq(2)
        
        # Enable tracking and delete file
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom2.txt")
        expect(last_response.status).to eq(200)
        
        # Should still correctly identify as orphaned (duplicate links don't break logic)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned.length).to eq(0)
      end

      it 'backfiller creates unlink event, then ledger creates link event for same datablock' do
        # Scenario: Orphaned datablock gets backfilled as unlinked, then new file reuses it
        disable_all_ledger_events
        
        # Create and delete file without tracking
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        original_datablock_id = wisdom_file.data_block_id
        
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Run backfill to create unlink event for orphaned datablock
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(orphaned: true)
        
        # Enable tracking and upload new file with same content (reuses datablock)
        enable_all_ledger_events
        new_wisdom_file = upload_file_via_api('athena', 'wisdom2.txt', WISDOM)
        expect(new_wisdom_file.data_block_id).to eq(original_datablock_id)
        
        # Verify new link event was created
        new_link_event = Metis::DataBlockLedger.where(
          file_id: new_wisdom_file.id,
          event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
        ).first
        expect(new_link_event).to be_present
        
        # Should NOT be orphaned (file currently exists with this datablock)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned).to be_empty
      end

      it 'backfiller creates unlink event, then ledger creates link and unlink events' do
        # Scenario: Orphaned datablock gets backfilled as unlinked, then new file reuses it and gets deleted
        disable_all_ledger_events
        
        # Create and delete file without tracking
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        original_datablock_id = wisdom_datablock.id
        
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Run backfill to create unlink event for orphaned datablock
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(orphaned: true)
        
        # Enable tracking and upload new file with same content (reuses datablock, creates link event)
        enable_all_ledger_events
        wisdom2_file = upload_file_via_api('athena', 'wisdom2.txt', WISDOM)
        expect(wisdom2_file.data_block_id).to eq(original_datablock_id)
        
        # Delete the new file (creates ledger unlink event)
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom2.txt")
        expect(last_response.status).to eq(200)
        
        # Should correctly identify as orphaned (linked and unlinked, no file exists)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned.length).to eq(1)
        expect(orphaned.first.id).to eq(wisdom_datablock.id)
      end

      it 'ledger creates link event, then backfiller runs and skips duplicate' do
        # Scenario: Tracking is on, file uploaded, then backfill runs with links: true
        enable_all_ledger_events
        
        # Upload file with tracking enabled
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        # Verify ledger link event was created
        ledger_link = Metis::DataBlockLedger.where(
          file_id: wisdom_file.id,
          event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
        ).exclude(triggered_by: Metis::DataBlockLedger::SYSTEM_BACKFILL).first
        expect(ledger_link).to be_present
        
        # Run backfill (should skip this file)
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(project_name: 'athena', links: true)
        
        # Should only have 1 link event (backfiller skipped)
        link_count = Metis::DataBlockLedger.where(
          file_id: wisdom_file.id,
          event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
        ).count
        expect(link_count).to eq(1)
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned).to be_empty
      end

      it 'ledger creates unlink event, then backfiller runs and skips duplicate' do
        # Scenario: Tracking is on, file deleted, then backfill runs with orphaned: true
        enable_all_ledger_events
        
        # Upload and delete file with tracking enabled
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        wisdom_datablock = wisdom_file.data_block
        
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
        
        # Verify ledger unlink event was created
        ledger_unlink = Metis::DataBlockLedger.where(
          data_block_id: wisdom_datablock.id,
          event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
        ).exclude(triggered_by: Metis::DataBlockLedger::SYSTEM_BACKFILL).first
        expect(ledger_unlink).to be_present
        
        # Run backfill (should skip this orphaned datablock)
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(orphaned: true)
        
        # Should only have 1 unlink event (backfiller skipped)
        unlink_count = Metis::DataBlockLedger.where(
          data_block_id: wisdom_datablock.id,
          event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
        ).count
        expect(unlink_count).to eq(1)
        
        # Should still correctly identify as orphaned
        orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
        expect(orphaned.length).to eq(1)
      end
    end
  end
end
