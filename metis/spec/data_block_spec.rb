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

  context '#check' do
    it 'checks for the existence of data blocks by md5' do
      wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
      stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

      wisdom_md5 = Digest::MD5.hexdigest(WISDOM)
      folly_md5 = Digest::MD5.hexdigest(WISDOM.reverse)

      token_header(:viewer)
      json_post('/check', md5s: [
        wisdom_md5,
        folly_md5
      ])

      expect(last_response.status).to eq(200)
      expect(json_body).to eq(found: [ wisdom_md5 ], missing: [ folly_md5])
    end
  end
end
