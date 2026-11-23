require_relative '../lib/remote_manager'

RSpec.describe Vulcan::RemoteManager do
  let(:ssh_pool) { double }
  let(:manager) { Vulcan::RemoteManager.new(ssh_pool) }

  describe '#read_file_to_memory' do
    it 'reads a remote file into memory' do
    end
  end
end
