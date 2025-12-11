require_relative '../lib/commands'

describe 'Metis Commands' do
  let(:metis_instance) { double('Metis') }

  describe Metis::Migrate do
    subject(:migrate) { described_class.new.execute(version: version) }
    let(:directory) { "db/migrations" }

    before do
      Sequel.extension(:migration)
      allow(Metis.instance).to receive(:db).and_return(metis_instance)
    end

    describe 'when a version is specified' do
      let(:version) { '001' }

      before do
        allow(Sequel::Migrator).to receive(:run).with(metis_instance, directory, target: version.to_i).and_return(true)
      end

      it 'calls run with a version number' do
        migrate

        expect(Sequel::Migrator)
        .to have_received(:run)
        .with(metis_instance, directory, target: version.to_i)
        .once
      end
    end

    describe 'without a version specified' do
      let(:version) { nil }

      before do
        allow(Sequel::Migrator).to receive(:run).with(metis_instance, directory).and_return(true)
      end

      it 'calls run without a version number' do
        migrate

        expect(Sequel::Migrator)
        .to have_received(:run)
        .with(metis_instance, directory)
        .once
      end
    end
  end

  describe Metis::Console do
    subject(:console) { described_class.new.execute }
    before do
      require 'irb'
      allow(ARGV).to receive(:clear)
      allow(IRB).to receive(:start)
    end

    it 'calls ARGV and IRB' do
      console

      expect(ARGV).to have_received(:clear).once
      expect(IRB).to have_received(:start).once
    end
  end

  describe Metis::CreateDb do
    let(:createdb_instance) { described_class.new }
    let(:config_double) { double('config') }
    let(:expected) { "Database is setup. Please run `bin/metis migrate `.\n" }
    subject(:create_db) { createdb_instance.execute() }

    before do
      allow(Metis).to receive(:instance).and_return(metis_instance)
      allow(metis_instance).to receive(:config).with(:db).and_return({ database: 'database' })
    end

    describe 'with @no_db = true' do
      before do
        allow(metis_instance).to receive(:setup_db).and_raise Sequel::DatabaseConnectionError
        allow(metis_instance).to receive(:configure).with(config_double)
        createdb_instance.setup(config_double)
        allow(createdb_instance).to receive(:create_db)
      end

      it 'calls create_db' do
        expect {
        create_db
        }.to output(expected).to_stdout

        expect(createdb_instance).to have_received(:create_db).once
      end
    end
  end
end

describe Metis::BackfillDataBlockLedger do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    default_bucket('athena')
    default_bucket('labors')
    default_bucket('backup')

    @metis_uid = Metis.instance.sign.uid

    set_cookie "#{Metis.instance.config(:metis_uid_name)}=#{@metis_uid}"

    @user = Etna::User.new(AUTH_USERS[:editor])
  end

  after(:each) do
    stubs.clear
  end

  context 'backfill ledger link' do
    before(:each) do
      token_header(:editor)
    end

    it "creates link events for all existing files" do
      disable_all_ledger_events
      
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      helmet_file = upload_file_via_api('athena', 'helmet.jpg', HELMET)

      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
      Metis::BackfillDataBlockLedger.new.execute(project_name: 'athena', links: true)

      # Should have 2 link_file_to_datablock events (one per file)
      link_events = Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).all
      expect(link_events.count).to eq(2)
      
      # Verify wisdom file ledger entry
      wisdom_ledger = Metis::DataBlockLedger.where(file_id: wisdom_file.id, event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).first
      expect(wisdom_ledger).to be_present
      expect(wisdom_ledger.project_name).to eq('athena')
      expect(wisdom_ledger.md5_hash).to eq(wisdom_file.data_block.md5_hash)
      expect(wisdom_ledger.triggered_by).to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
      
      # Verify helmet file ledger entry
      helmet_ledger = Metis::DataBlockLedger.where(file_id: helmet_file.id, event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).first
      expect(helmet_ledger).to be_present
      expect(helmet_ledger.project_name).to eq('athena')
      expect(helmet_ledger.md5_hash).to eq(helmet_file.data_block.md5_hash)
      expect(helmet_ledger.triggered_by).to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
      
      # Clean up
      wisdom_file.delete
      helmet_file.delete
    end

    it "skips files already in ledger when backfilling links" do
      disable_all_ledger_events
      
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)

      backfill_ledger = Metis::BackfillDataBlockLedger.new
      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')

      # Run backfill first time - should create link event
      backfill_ledger.execute(project_name: 'athena', links: true)
      expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).count).to eq(1)

      # Run backfill second time - should NOT create duplicate link event
      backfill_ledger.execute(project_name: 'athena', links: true)

      # Should still be 1 link_file_to_datablock event, not duplicated
      expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).count).to eq(1)
      
      wisdom_file.delete
    end

    it "creates links events across all projects" do
      disable_all_ledger_events
      
      # Create files with the same data block
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
      backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)

      backfill_ledger = Metis::BackfillDataBlockLedger.new
      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
      backfill_ledger.execute(project_name: 'athena', links: true)
      backfill_ledger.execute(project_name: 'labors', links: true)
      backfill_ledger.execute(project_name: 'backup', links: true)
      expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).count).to eq(3)
    end

    it 'skips files that already have tracked link events' do
      # Create a untracked file
      disable_all_ledger_events

      untracked_file = upload_file_via_api('athena', 'untracked1.txt', WISDOM)

      # Enable ledger tracking (allow real ledger events)
      enable_all_ledger_events
    
      # Create a tracked file
      tracked_file = upload_file_via_api('athena', 'tracked2.txt', HELMET)
      
      # Run backfill - should skip tracked files, but backfill the others
      backfill_ledger = Metis::BackfillDataBlockLedger.new
      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
      backfill_ledger.execute(project_name: 'athena', links: true)
      
      # Verify tracked_file1 (created before ledger enabled) now has backfilled link event
      tracked_files = Metis::DataBlockLedger.where(
        event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK,
        triggered_by: Metis::DataBlockLedger::SYSTEM_BACKFILL,
      ).all
      expect(tracked_files.count).to eq(1)
      expect(tracked_files.first.triggered_by).to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
    end
  end

  context 'backfill ledger unlink' do

    # We typical run the link command first, so this is why it is being invoked here.

    it "creates unlink events for orphaned datablocks" do
      disable_all_ledger_events
      
      # Create two files
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      helmet_file = upload_file_via_api('athena', 'helmet.jpg', HELMET)
      
      wisdom_datablock_id = wisdom_file.data_block_id
      helmet_datablock_id = helmet_file.data_block_id

      # Delete both files using the API (simulating real user deletion)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)
      
      delete("/athena/file/remove/files/helmet.jpg")
      expect(last_response.status).to eq(200)


      # Mock the ask_user method to return 'y' automatically
      backfill_ledger = Metis::BackfillDataBlockLedger.new

      # Run backfill should detect both orphaned datablocks
      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
      backfill_ledger.execute(project_name: 'athena', links: true)
      backfill_ledger.execute(orphaned: true)

      # There should be no link events for the datablocks
      link_events = Metis::DataBlockLedger.where(
        event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
      ).all
      expect(link_events.count).to eq(0)
      
      # Should have created unlink events for both datablocks
      unlink_events = Metis::DataBlockLedger.where(
        event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
      ).all

      expect(unlink_events.count).to eq(2)
      
      # Verify wisdom file unlink event
      wisdom_unlink = unlink_events.find { |e| e.data_block_id == wisdom_datablock_id }
      expect(wisdom_unlink).to be_present
      expect(wisdom_unlink.project_name).to be_nil
      expect(wisdom_unlink.triggered_by).to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
      expect(wisdom_unlink.size).to eq(WISDOM.length)
      
      # Verify helmet file unlink event
      helmet_unlink = unlink_events.find { |e| e.data_block_id == helmet_datablock_id }
      expect(helmet_unlink).to be_present
      expect(helmet_unlink.project_name).to be_nil
      expect(helmet_unlink.triggered_by).to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
      expect(helmet_unlink.size).to eq(HELMET.length)
    end

    it "skips files already in ledger when backfilling orphaned datablocks" do
      disable_all_ledger_events
      
      wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
      
      wisdom_datablock_id = wisdom_file.data_block_id

      # Delete the file using the API (simulating real user deletion)
      token_header(:editor)
      delete("/athena/file/remove/files/wisdom.txt")
      expect(last_response.status).to eq(200)

      backfill_ledger = Metis::BackfillDataBlockLedger.new
      allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')

      # Run backfill first time - should create unlink event
      backfill_ledger.execute(project_name: 'athena', links: true)
      backfill_ledger.execute(orphaned: true)
      expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).count).to eq(0)
      expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK).count).to eq(1)

      # Run backfill second time - should NOT create duplicate unlink event
      backfill_ledger.execute(orphaned: true)

      # Should still be 1 unlink_file_from_datablock event, not duplicated
      expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK).count).to eq(1)
    end

    it "skips datablocks when they are linked to multiple projects" do
        disable_all_ledger_events
        
        # Create files with the same data block
        wisdom_file = upload_file_via_api('athena', 'wisdom.txt', WISDOM)
        labors_file = upload_file_via_api('labors', 'labors.txt', WISDOM)
        backup_file = upload_file_via_api('backup', 'backup.txt', WISDOM)
  
        token_header(:editor)
        delete("/athena/file/remove/files/wisdom.txt")
        expect(last_response.status).to eq(200)
  
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
  
        # Run backfill should detect no orphaned datablocks
        backfill_ledger.execute(project_name: 'athena', links: true)
        backfill_ledger.execute(project_name: 'labors', links: true)
        backfill_ledger.execute(project_name: 'backup', links: true)
        backfill_ledger.execute(orphaned: true)
  
        expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK).count).to eq(2)
        expect(Metis::DataBlockLedger.where(event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK).count).to eq(0)
    end

    it 'skips files that already have tracked link events' do
        disable_all_ledger_events
  
        # Create a untracked file and delete it
        untracked_file = upload_file_via_api('athena', 'untracked1.txt', WISDOM)
        token_header(:editor)
        delete("/athena/file/remove/files/untracked1.txt")
  
        # Enable ledger tracking (allow real ledger events)
        enable_all_ledger_events
      
        # Create a tracked file
        tracked_file = upload_file_via_api('athena', 'tracked2.txt', HELMET)
        delete("/athena/file/remove/files/tracked2.txt")
        
        # Run backfill - should skip tracked files, but backfill the others
        backfill_ledger = Metis::BackfillDataBlockLedger.new
        allow_any_instance_of(Metis::BackfillDataBlockLedger).to receive(:ask_user).and_return('y')
        backfill_ledger.execute(orphaned: true)
        
        # Verify tracked_file1 (created before ledger enabled) now has backfilled link event
        tracked_files = Metis::DataBlockLedger.where(
          triggered_by: Metis::DataBlockLedger::SYSTEM_BACKFILL,
          event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK,
        ).all
        expect(tracked_files.count).to eq(1)
        expect(tracked_files.first.triggered_by).to eq(Metis::DataBlockLedger::SYSTEM_BACKFILL)
      end
    end
end