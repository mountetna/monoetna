require_relative '../lib/commands'

describe 'Vesta Commands' do
  let(:vesta_instance) { double('Vesta') }
  let(:vesta_db_instance) { double('Sequel::Database') }

  describe Vesta::CollectDLStats do
    subject(:collect_global_stats) {
      described_class.new
    }

    it "ddsa" do
    end
  end

  describe Vesta::Migrate do
    subject(:migrate) { described_class.new.execute(version: version) }
    let(:directory) { "db/migrations" }

    before do
      Sequel.extension(:migration)
      allow(Vesta).to receive(:instance).and_return(vesta_instance)
      allow(Vesta.instance).to receive(:db).and_return(vesta_db_instance)
    end

    describe 'when a version is specified' do
      let(:version) { '001' }

      before do
        allow(Sequel::Migrator).to receive(:run).with(vesta_db_instance, directory, target: version.to_i).and_return(true)
      end

      it 'calls run with a version number' do
        migrate

        expect(Sequel::Migrator)
        .to have_received(:run)
        .with(vesta_db_instance, directory, target: version.to_i)
        .once
      end
    end

    describe 'without a version specified' do
      let(:version) { nil }

      before do
        allow(Sequel::Migrator).to receive(:run).with(vesta_db_instance, directory).and_return(true)
      end

      it 'calls run without a version number' do
        migrate

        expect(Sequel::Migrator)
        .to have_received(:run)
        .with(vesta_db_instance, directory)
        .once
      end
    end
  end

  describe Vesta::Console do
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
end
