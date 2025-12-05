require 'digest'
require_relative '../lib/commands'

describe Metis::DataBlock do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    default_bucket('athena')

    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"

    stub_event_log(:athena)
  end

  after(:each) do
    stubs.clear
  end


  context '#compute_hash!' do
    def validate_readable_block(block)
      expect(block.has_data?).to be_truthy

      unless block.temp_hash?
        expect(Metis::File.md5(block.location)).to eql(block.md5_hash)
      end
    end

    def validate_block_converged(block)
      return if Metis::DataBlock.where(id: block.id).first.nil?
      return if block.removed
      return unless block.temp_hash?

      validate_readable_block(block)
    end

    def self.test_concurrent_action(&block)
      before(:each) do
        new_db.compute_hash!(&block)
      end

      it "converges on an end state" do
        validate_block_converged(new_db)
        validate_readable_block(new_file.reload.data_block)
        expect(new_db.temp_hash?).to be_falsey
      end
    end

    def self.test_with_processing_error(&block)
      before(:each) do
        expect do
          new_db.compute_hash!(&block)
        end.to raise_error("Test")

        expect(new_db.temp_hash?).to be_truthy
      end

      it 'preserves a working file after the failure' do
        validate_readable_block(new_file.reload.data_block)
      end

      it "converges on its next run on a valid file" do
        new_db.compute_hash!
        validate_block_converged(new_db)
        validate_readable_block(new_file.reload.data_block)
        expect(new_file.data_block.temp_hash?).to be_falsey
        expect(::File.symlink?(new_db.location)).to be_falsey
      end
    end

    describe "processing a new temp file" do

      let(:temp_hash) { "#{Metis::DataBlock::TEMP_PREFIX}abcdefg" }
      let(:new_file) do
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM, temp_hash)
        create_file('athena', 'wisdom.txt', WISDOM, md5_hash: temp_hash)
      end

      let(:expected_hash) do
        Digest::MD5.hexdigest(WISDOM)
      end

      let(:new_db) do
        new_file.data_block
      end

      let(:existing_file) do
        stubs.create_file('athena', 'files', 'wisdom2.txt', WISDOM, expected_hash)
        create_file('athena', 'wisdom2.txt', WISDOM, md5_hash: expected_hash)
      end

      describe 'for a new hash' do
        describe 'when a data block by the same hash concurrently is created' do
          test_concurrent_action do |test_cmd|
            if test_cmd.first == :sym
              # create the existing file, at this very moment.
              ::File.write(test_cmd[2], WISDOM)
            end
          end
        end

        describe 'when the update to the block md5 fails after the sym link' do
          test_with_processing_error do |test_cmd|
            if test_cmd.first == :new_update
              raise "Test"
            end
          end
        end
      end

      describe "with existing block pointing to same hash" do
        before(:each) { existing_file }

        it 'converges sanely' do
          new_db.compute_hash!

          validate_block_converged(new_db)
          validate_readable_block(existing_file.reload.data_block)
          validate_readable_block(new_file.reload.data_block)
        end

        it 'logs a REUSE_DATABLOCK ledger entry when deduplicating' do
          enable_all_ledger_events
          
          new_db.compute_hash!

          validate_block_converged(new_db)
          validate_readable_block(existing_file.reload.data_block)
          validate_readable_block(new_file.reload.data_block)

          # Assert log_deduplicate was called when compute_hash! found existing block
          dedupe_event = Metis::DataBlockLedger.where(
            event_type: Metis::DataBlockLedger::REUSE_DATABLOCK,
            data_block_id: existing_file.data_block_id
          ).first
          expect(dedupe_event).to be_present
          expect(dedupe_event.project_name).to eq('athena')
          expect(dedupe_event.triggered_by).to eq('system')
        end

        describe 'if something happens after temp delete but before database update' do
          test_with_processing_error do |test_cmd|
            if test_cmd.first == :temp_destroy
              raise "Test"
            end
          end

          it "The existing file is still ok" do
            validate_readable_block(existing_file.reload.data_block)
          end
        end

        describe 'if something happens before existing temp file deletion' do
          test_with_processing_error do |test_cmd|
            if test_cmd.first == :temp_delete
              raise "Test"
            end
          end

          it "The existing file is still ok" do
            validate_readable_block(existing_file.reload.data_block)
          end
        end
      end
    end

    it 'computes the md5 sum of the block if it is temporary (not yet computed)' do
      # We create the data block with a temporary hash assigned
      temp_hash = "temp-ef15c9bd4c7836612b1567f4c8396726"
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, md5_hash: temp_hash)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM, temp_hash)

      wisdom_data = wisdom_file.data_block

      # We expect some actual hash computation to happen
      expect(Metis::File).to receive(:md5).with(wisdom_data.location).and_call_original

      wisdom_data.compute_hash!

      wisdom_data.refresh

      # the hash is now the actual md5
      expect(wisdom_data.md5_hash).to eq(Digest::MD5.hexdigest(WISDOM))
      expect(wisdom_data.has_data?).to be_truthy
    end

    it 'does not recompute hashes' do
      # We initialize the data block with its actual hash
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      wisdom_data = wisdom_file.data_block

      expect(Metis::File).not_to receive(:md5)

      wisdom_data.compute_hash!

      wisdom_data.refresh
      expect(wisdom_data.md5_hash).to eq(Digest::MD5.hexdigest(WISDOM))
      expect(wisdom_data.has_data?).to be_truthy
    end
  end

  context '#backup!' do
    before(:each) do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      @wisdom_data = @wisdom_file.data_block
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
    end

    it 'backs up the file to AWS Glacier' do
      glacier_stub('metis-test-athena')
      @wisdom_data.backup!

      expect(Metis::DataBlock.count).to eq(1)

      @wisdom_data.refresh
      expect(@wisdom_data.archive_id).to eq('archive_id')
    end
  end

  context '#remove!' do
    it 'removes the data_block from disk and sets removed flag' do
      @creation_time = DateTime.now - 10
      Timecop.freeze(@creation_time)
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      wisdom_data = wisdom_file.data_block

      expect(::File.exists?(wisdom_data.location)).to eq(true)
      expect(wisdom_data.removed).to eq(false)
      expect(wisdom_data.updated_at.iso8601).to eq(@creation_time.to_s)

      @update_time = DateTime.now
      Timecop.freeze(@update_time)

      wisdom_data.remove!

      expect(::File.exists?(wisdom_data.location)).to eq(false)
      expect(wisdom_data.removed).to eq(true)
      expect(wisdom_data.updated_at.iso8601).to eq(@update_time.to_s)
      Timecop.return
    end

    it 'only changes removed flag if the block location does not exist on disk' do
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      wisdom_data = wisdom_file.data_block

      ::File.delete(wisdom_data.location)

      expect(::File.exists?(wisdom_data.location)).to eq(false)
      expect(wisdom_data.removed).to eq(false)

      wisdom_data.remove!

      expect(::File.exists?(wisdom_data.location)).to eq(false)
      expect(wisdom_data.removed).to eq(true)
    end

    it 'does not execute anything if block already removed' do
      past_time = DateTime.now - 10
      Timecop.freeze(past_time)
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      wisdom_data = wisdom_file.data_block
      wisdom_data.update(removed: true)

      expect(::File.exists?(wisdom_data.location)).to eq(true)
      expect(wisdom_data.removed).to eq(true)
      expect(wisdom_data.updated_at.iso8601).to eq(past_time.to_s)

      Timecop.return

      wisdom_data.remove!

      # Since no action should have been taken
      expect(::File.exists?(wisdom_data.location)).to eq(true)
      expect(wisdom_data.removed).to eq(true)
      expect(wisdom_data.updated_at.iso8601).to eq(past_time.to_s)
    end
  end
end

describe DataBlockController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    default_bucket('athena')

    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"
  end

  after(:each) do
    stubs.clear
  end

  context '#exists' do
    it 'checks for the existence of data blocks by md5' do
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
      folly_file = create_file('ate', 'folly.txt', WISDOM.reverse)
      stubs.create_file('ate', 'files', 'folly.txt', WISDOM.reverse)

      wisdom_md5 = Digest::MD5.hexdigest(WISDOM)
      folly_md5 = Digest::MD5.hexdigest(WISDOM.reverse)
      helmet_md5 = Digest::MD5.hexdigest(HELMET)

      token_header(:viewer)
      json_post('/api/exists', md5s: [
        wisdom_md5,
        folly_md5,
        helmet_md5
      ])

      expect(last_response.status).to eq(200)
      expect(json_body).to eq(found: [wisdom_md5, folly_md5], missing: [helmet_md5])
    end

    it 'finds data blocks only for a specified project' do
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)
      folly_file = create_file('ate', 'folly.txt', WISDOM.reverse)
      stubs.create_file('ate', 'files', 'folly.txt', WISDOM.reverse)

      wisdom_md5 = Digest::MD5.hexdigest(WISDOM)
      folly_md5 = Digest::MD5.hexdigest(WISDOM.reverse)

      token_header(:viewer)
      json_post('/api/exists', md5s: [
        wisdom_md5,
        folly_md5
      ], project_name: 'athena')

      expect(last_response.status).to eq(200)
      expect(json_body).to eq(found: [wisdom_md5], missing: [folly_md5])
    end
  end

  context '#vacuum_datablocks with existing backfilled records' do
    before(:each) do
      # Mock ledger methods during setup to prevent events from being created
    end

    it 'vacuums orphaned datablocks for backfilled records' do
      disable_all_ledger_events
      
      # Create and then delete files to make datablocks orphaned, then backfill
      result = athena_backfilled_lifecycle
      wisdom_data_block = Metis::DataBlock.where(id: result[:wisdom_data_block_id]).first
      helmet_data_block = Metis::DataBlock.where(id: result[:helmet_data_block_id]).first
      
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/backfilled', {})
      
      expect(last_response.status).to eq(200)
      response = json_body

      # Summary stats
      expect(response[:vacuumed].length).to eq(2)
      expect(response[:summary][:space_freed]).to eq(wisdom_data_block.size + helmet_data_block.size)
      
      # Wisdom data block
      expect(response[:vacuumed].first[:md5_hash]).to eq(wisdom_data_block.md5_hash)
      expect(response[:vacuumed].first[:size]).to eq(wisdom_data_block.size)
      expect(response[:vacuumed].first[:location]).to eq(wisdom_data_block.location)
      
      # Helmet data block
      expect(response[:vacuumed].last[:md5_hash]).to eq(helmet_data_block.md5_hash)
      expect(response[:vacuumed].last[:size]).to eq(helmet_data_block.size)
      expect(response[:vacuumed].last[:location]).to eq(helmet_data_block.location)

      # Verify datablock is marked as removed
      wisdom_data_block.reload
      helmet_data_block.reload
      expect(wisdom_data_block.removed).to be_truthy
      expect(helmet_data_block.removed).to be_truthy

      # Verify the file contents are no longer present
      expect(::File.exists?(wisdom_data_block.location)).to eq(false)
      expect(::File.exists?(helmet_data_block.location)).to eq(false)
      
      # Verify vacuum event was logged
      vacuum_event = Metis::DataBlockLedger.where(
        md5_hash: wisdom_data_block.md5_hash,
        event_type: Metis::DataBlockLedger::REMOVE_DATABLOCK
      ).first
      expect(vacuum_event).to be_present

      vacuum_event = Metis::DataBlockLedger.where(
        md5_hash: helmet_data_block.md5_hash,
        event_type: Metis::DataBlockLedger::REMOVE_DATABLOCK
      ).first
      expect(vacuum_event).to be_present
    end

  end

  context '#vacuum_datablocks' do
    it 'vacuums orphaned datablocks for a project' do
      enable_all_ledger_events
      
      # Create file via upload API (triggers log_link automatically)
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      
      # Delete the file via API (triggers log_unlink automatically)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      
      # Verify unlink event was created (not SYSTEM_BACKFILL)
      unlink_event = Metis::DataBlockLedger.where(
        data_block_id: wisdom_file.data_block_id,
        event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
      ).first
      expect(unlink_event).to be_present
      expect(unlink_event.triggered_by).not_to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
      
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/athena', {})
      
      expect(last_response.status).to eq(200)
      response = json_body
      expect(response[:vacuumed].length).to eq(1)
      expect(response[:vacuumed].first[:md5_hash]).to eq(wisdom_file.data_block.md5_hash)
      expect(response[:errors]).to be_empty
      expect(response[:summary][:total_vacuumed]).to eq(1)
      expect(response[:summary][:space_freed]).to be > 0
      
      # Verify datablock is marked as removed
      wisdom_file.data_block.reload
      expect(wisdom_file.data_block.removed).to be_truthy
      
      # Assert log_vacuum was called and vacuum event was logged
      vacuum_event = Metis::DataBlockLedger.where(
        md5_hash: wisdom_file.data_block.md5_hash,
        event_type: Metis::DataBlockLedger::REMOVE_DATABLOCK
      ).first
      expect(vacuum_event).to be_present
      expect(vacuum_event.project_name).to eq('athena')
      expect(vacuum_event.data_block_id).to eq(wisdom_file.data_block_id)
    end

    it 'only vacuums datablocks for the specified project' do
      # Create files in two projects via upload API
      default_bucket('ate')
      athena_file = upload_file_via_api('athena', 'athena.txt', WISDOM)
      ate_file = upload_file_via_api('ate', 'ate.txt', HELMET)
      
      # Delete athena file via API (triggers log_unlink automatically)
      token_header(:editor)
      delete("/athena/file/remove/files/athena.txt")
      expect(last_response.status).to eq(200)
      
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/athena', {})
      
      expect(last_response.status).to eq(200)
      response = json_body
      
      # Should only vacuum athena's orphaned blocks, not ate's
      expect(response[:vacuumed].length).to eq(1)
      expect(response[:summary][:project_name]).to eq('athena')
      
      # Assert log_vacuum was called for athena's datablock
      vacuum_event = Metis::DataBlockLedger.where(
        data_block_id: athena_file.data_block_id,
        event_type: Metis::DataBlockLedger::REMOVE_DATABLOCK
      ).first
      expect(vacuum_event).to be_present
      expect(vacuum_event.project_name).to eq('athena')
    end

    it 'does not vacuum datablocks still referenced by files' do
      # Create file via upload API (triggers log_link automatically)
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/athena', {})
      
      expect(last_response.status).to eq(200)
      response = json_body
      expect(response[:vacuumed]).to be_empty
      expect(response[:summary][:total_vacuumed]).to eq(0)
    end

    it 'does not vacuum already vacuumed datablocks' do
      # Create file via upload API (triggers log_link automatically)
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      
      # Delete file via API (triggers log_unlink automatically)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      
      # Vacuum once
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/athena', {})
      expect(last_response.status).to eq(200)
      
      # Try to vacuum again - should find nothing
      json_post('/api/vacuum_datablocks/athena', {})
      expect(last_response.status).to eq(200)
      response = json_body
      expect(response[:vacuumed]).to be_empty
    end

    it 'handles multiple orphaned datablocks' do
      # Create files via upload API (triggers log_link automatically)
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      helmet_file = upload_file_via_api('athena', 'helmet.jpg', HELMET)
      
      # Delete files via API (triggers log_unlink automatically)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      delete("/athena/file/remove/files/helmet.jpg")
      expect(last_response.status).to eq(200)
      
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/athena', {})
      
      expect(last_response.status).to eq(200)
      response = json_body
      expect(response[:vacuumed].length).to eq(2)
      expect(response[:summary][:total_vacuumed]).to eq(2)
    end

    it 'requires admin permissions' do
      token_header(:editor)
      json_post('/api/vacuum_datablocks/athena', {})
      
      expect(last_response.status).to eq(403)
    end

    it 'returns empty result when no orphaned datablocks exist' do
      token_header(:supereditor)
      json_post('/api/vacuum_datablocks/athena', {})
      
      expect(last_response.status).to eq(200)
      response = json_body
      expect(response[:vacuumed]).to be_empty
      expect(response[:summary][:total_vacuumed]).to eq(0)
      expect(response[:summary][:space_freed]).to eq(0)
    end
  end

end
