require 'digest'

describe Metis::DataBlock do
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

  context '#compute_hash!' do
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
end
