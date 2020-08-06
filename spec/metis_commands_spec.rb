describe 'Metis Commands' do
  describe Metis::DeleteOrphanDataBlocks do
    subject(:delete_orphan_data_blocks) { described_class.new.execute }

    it "does not delete used data blocks" do
      expected = "Found 0 orphaned data blocks to be deleted.\n"

      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
      stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)

      expect {
        delete_orphan_data_blocks
      }.to output(expected).to_stdout

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)

      # Clean up the test
      @wisdom_file.delete
      @helmet_file.delete
    end

    it "deletes orphaned data blocks" do
      @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
      stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

      expected = "Found 1 orphaned data blocks to be deleted.\nDeleted data_block with hash #{@wisdom_file.data_block.md5_hash}\n"

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(2)

      @wisdom_file.update({data_block: @helmet_file.data_block})

      expect {
        delete_orphan_data_blocks
      }.to output(expected).to_stdout

      expect(Metis::File.count).to eq(2)
      expect(Metis::DataBlock.count).to eq(1)

      # Clean up the test
      @wisdom_file.delete
      @helmet_file.delete
    end

    it "does not delete the zero-hash data block" do
        expected = "Found 0 orphaned data blocks to be deleted.\n"

        zero_hash = 'd41d8cd98f00b204e9800998ecf8427e'

        @wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
        stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

        @helmet_file = create_file('athena', 'helmet.jpg', HELMET)
        stubs.create_file('athena', 'files', 'helmet.jpg', HELMET)

        @zero_hash_data_block = create(:data_block,
          description: 'zero-byte hash',
          md5_hash: zero_hash
        )

        expect(Metis::File.count).to eq(2)
        expect(Metis::DataBlock.count).to eq(3)

        expect {
          delete_orphan_data_blocks
        }.to output(expected).to_stdout

        expect(Metis::File.count).to eq(2)
        expect(Metis::DataBlock.count).to eq(3)

        # Clean up the test
        @wisdom_file.delete
        @helmet_file.delete
        @zero_hash_data_block.delete
      end
  end

  let(:metis_instance) { double('Metis') }

  describe Metis::Help do
    subject(:help) { described_class.new.execute }

    let(:command_double) { double('command', usage: 'Usage') }
    let(:expected) { "Commands:\nUsage\n" }

    before do
      allow(Metis).to receive(:instance).and_return(metis_instance)
      allow(metis_instance).to receive(:commands).and_return({"test" => command_double})
    end

    it "calls puts once for each command present" do
      expect {
        help
      }.to output(expected).to_stdout
    end
  end

  describe Metis::Migrate do
    subject(:migrate) { described_class.new.execute(version) }
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
