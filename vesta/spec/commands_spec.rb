require_relative '../lib/commands'

describe 'Vesta Commands' do
  describe Vesta::CollectDLStats do
    subject(:collect_global_stats) {
      described_class.new
    }

    it "collects data from other apps" do
      command = Vesta::CollectDLStats.new

      command.execute
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
