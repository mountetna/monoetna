require_relative '../lib/commands'

describe 'Metis Commands' do
  describe Metis::RemoveOrphanDataBlocks do
    subject(:remove_orphan_data_blocks) { described_class.new.execute }

    after(:each) do
      stubs.clear
    end

    it "does not remove used data blocks" do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
      stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)
      wisdom_file_block_location = @wisdom_file.data_block.location
      helmet_file_block_location = @helmet_file.data_block.location
      expect(::File.exists?(wisdom_file_block_location)).to eq(true)
      expect(::File.exists?(helmet_file_block_location)).to eq(true)

      remove_orphan_data_blocks

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)
      expect(::File.exists?(wisdom_file_block_location)).to eq(true)
      expect(::File.exists?(helmet_file_block_location)).to eq(true)

      # Clean up the test
      @wisdom_file.delete
      @helmet_file.delete
    end

    it "removes orphaned data blocks" do
      glacier_stub('metis-test-athena')

      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
      stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)

      wisdom_block_md5_hash = @wisdom_file.data_block.md5_hash
      wisdom_file_block_location = @wisdom_file.data_block.location
      helmet_file_block_location = @helmet_file.data_block.location
      @wisdom_file.data_block.update(archive_id: 'archive_id')

      expect(::File.exists?(wisdom_file_block_location)).to eq(true)
      expect(::File.exists?(helmet_file_block_location)).to eq(true)

      @wisdom_file.update({data_block: @helmet_file.data_block})

      remove_orphan_data_blocks

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)
      expect(::File.exists?(wisdom_file_block_location)).to eq(false)
      expect(::File.exists?(helmet_file_block_location)).to eq(true)

      removed_data_block = Metis::DataBlock.find({:md5_hash => wisdom_block_md5_hash})

      expect(removed_data_block.removed).to eq(true)

      # Clean up the test
      @wisdom_file.delete
      @helmet_file.delete
    end

    it "does not remove already removed data blocks" do
      past_time = DateTime.now - 10
      Timecop.freeze(past_time)
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
      stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)

      wisdom_block = @wisdom_file.data_block
      wisdom_file_block_location = @wisdom_file.data_block.location
      helmet_file_block_location = @helmet_file.data_block.location

      expect(::File.exists?(wisdom_file_block_location)).to eq(true)
      expect(::File.exists?(helmet_file_block_location)).to eq(true)

      wisdom_block.update(removed: true)
      @wisdom_file.data_block = wisdom_block

      Timecop.return

      remove_orphan_data_blocks

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)

      expect(wisdom_block.updated_at.iso8601).to eq(past_time.to_s)

      # Clean up the test
      @wisdom_file.delete
      @helmet_file.delete
    end

    it "does not remove the zero-hash data block" do
        zero_hash = 'd41d8cd98f00b204e9800998ecf8427e'

        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
        stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

        @zero_hash_data_block = create(:data_block,
          description: 'zero-byte hash',
          md5_hash: zero_hash,
          size: 0
        )
        stubs.create_file('athena', 'files', @zero_hash_data_block.md5_hash, '')

        wisdom_file_block_location = @wisdom_file.data_block.location
        helmet_file_block_location = @helmet_file.data_block.location
        zero_hash_block_location = @zero_hash_data_block.location

        expect(Metis::File.count).to eq(2)
        expect(Metis::DataBlock.count).to eq(3)

        expect(::File.exists?(wisdom_file_block_location)).to eq(true)
        expect(::File.exists?(helmet_file_block_location)).to eq(true)
        expect(::File.exists?(zero_hash_block_location)).to eq(true)

        remove_orphan_data_blocks

        expect(Metis::File.count).to eq(2)
        expect(Metis::DataBlock.count).to eq(3)

        expect(::File.exists?(wisdom_file_block_location)).to eq(true)
        expect(::File.exists?(helmet_file_block_location)).to eq(true)
        expect(::File.exists?(zero_hash_block_location)).to eq(true)

        # Clean up the test
        @wisdom_file.delete
        @helmet_file.delete
        @zero_hash_data_block.remove!
        @zero_hash_data_block.delete
      end
  end

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
