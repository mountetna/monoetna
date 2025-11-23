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
    
    # Enable ledger for these tests
    ENV['METIS_LEDGER_ENABLED'] = 'true'
  end

  after(:each) do
    stubs.clear
  end

  describe '.find_orphaned_datablocks' do
    it 'returns orphaned datablocks that were linked and then unlinked' do
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

    it 'returns empty array when no orphaned datablocks exist' do
      orphaned = Metis::DataBlockLedger.find_orphaned_datablocks('athena')
      expect(orphaned).to be_empty
    end

    describe 'cross project' do
      it 'does not return datablocks that are still in use by another project' do
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

    it 'handles datablocks with multiple link and unlink events correctly' do
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

  end

  describe 'find_orphaned_datablocks_legacy' do
    it 'returns orphaned datablocks' do
      # Disable ledger for backfill tests
      ENV['METIS_LEDGER_ENABLED'] = 'false'
      
      # Create files
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      helmet_file = create_file('athena', 'helmet.jpg', HELMET)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)
      
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
      backfill_ledger.execute('-orphaned')
      
      # Find orphaned legacy datablocks
      orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_legacy
      
      expect(orphaned.length).to eq(2)
      orphaned_ids = orphaned.map(&:id)
      expect(orphaned_ids).to include(wisdom_datablock.id, helmet_datablock.id)
    end

    it 'does not return datablocks that have been vacuumed' do
      # Disable ledger for backfill tests
      ENV['METIS_LEDGER_ENABLED'] = 'false'
      
      # Create a file
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
      wisdom_datablock = wisdom_file.data_block
      
      # Delete file via API
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      
      # Run backfill to create SYSTEM_BACKFILL unlink event
      backfill_ledger = Metis::BackfillDataBlockLedger.new
      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
      backfill_ledger.execute('-orphaned')
      
      # Vacuum the datablock via legacy API (marks it as removed and creates REMOVE_DATABLOCK event)
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/legacy', {})
      expect(last_response.status).to eq(200)
      
      # Should not be orphaned because it's been vacuumed (marked as removed)
      orphaned = Metis::DataBlockLedger.find_orphaned_datablocks_legacy
      
      expect(orphaned).to be_empty
    end
  end
end
