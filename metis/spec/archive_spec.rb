require 'digest/md5'
require_relative '../lib/commands'

describe Metis::ChecksumFiles do
  before(:each) do
    default_bucket('athena')
    @cmd = Metis::ChecksumFiles.new
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  it 'calculates md5s for new data blocks' do
    # a new data block has a temporary hash, recognizable by its pattern
    temp_hash = "temp-ef15c9bd4c7836612b1567f4c8396726"
    wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, md5_hash: temp_hash)
    stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM, temp_hash)
    wisdom_data = wisdom_file.data_block

    # skip archiving
    #wisdom_data.update(archive_id: 'skip_archiving')
    #wisdom_data.refresh

    @cmd.execute

    wisdom_data.refresh

    # now the hash is updated, and the block has been moved to its new location on disk
    expect(wisdom_data.md5_hash).to eq(Digest::MD5.hexdigest(WISDOM))
    expect(wisdom_data).to be_has_data
  end

  it "does not try to checksum removed data blocks" do
    wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
    stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

    temp_md5_hash = "temp-#{wisdom_file.data_block.md5_hash}"
    wisdom_file.data_block.update(md5_hash: temp_md5_hash, removed: true)

    expected = "Found 0 data blocks to be checksummed.\n"

    expect {
      @cmd.execute
    }.to output(expected).to_stdout

    expect(wisdom_file.data_block.md5_hash).to eq(temp_md5_hash)
  end

  it 'consolidates duplicate MD5s' do
    # a newly created file which needs to be hashed
    temp_hash = "temp-ef15c9bd4c7836612b1567f4c8396726"
    wisdom_file = create_file('athena', 'wisdom.txt', WISDOM, md5_hash: temp_hash)
    stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM, temp_hash)
    wisdom_data = wisdom_file.data_block

    # an existing file with the same md5
    old_wisdom_file = create_file('athena', 'old_wisdom.txt', WISDOM)
    stubs.create_file('athena', 'files', 'old_wisdom.txt', WISDOM)
    old_wisdom_data = old_wisdom_file.data_block

    # an innocent bystander file
    helmet_file = create_file('athena', 'helmet.txt', HELMET)
    stubs.create_file('athena', 'files', 'helmet.txt', HELMET)
    helmet_data = helmet_file.data_block

    # skip archiving
    old_wisdom_data.update(archive_id: 'skip_archiving')
    wisdom_data.update(archive_id: 'skip_archiving')
    helmet_data.update(archive_id: 'skip_archiving')

    @cmd.execute

    # there are only two data blocks; they have data
    expect(Metis::DataBlock.count).to eq(2)
    expect(Metis::DataBlock.all).to all(be_has_data)

    # the duplicate is gone, along with its data
    expect(Metis::DataBlock[wisdom_data.id]).to be_nil
    expect(wisdom_data).not_to be_has_data

    [ old_wisdom_file, old_wisdom_data,
      wisdom_file, helmet_file, helmet_data ].each(&:refresh)

    # there are still three files, pointing to the correct blocks
    expect(Metis::File.count).to eq(3)
    expect(old_wisdom_file.data_block_id).to eq(old_wisdom_data.id)
    expect(wisdom_file.data_block_id).to eq(old_wisdom_data.id)
    expect(helmet_file.data_block_id).to eq(helmet_data.id)

  end
end

describe Metis::Archive do
  before(:each) do
    default_bucket('athena')
    @cmd = Metis::Archive.new
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end


  it 'archives un-archived data blocks' do
    glacier_stub('metis-test-athena')
    wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
    stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

    wisdom_data = wisdom_file.data_block

    @cmd.execute

    wisdom_data.refresh
    expect(wisdom_data.archive_id).to eq('archive_id')
  end

  it "does not try to archive removed data blocks" do
    wisdom_file = create_file('athena', 'wisdom.txt', WISDOM)
    stubs.create_file('athena', 'files', 'wisdom.txt', WISDOM)

    expect(wisdom_file.data_block.archive_id).to eq(nil)
    wisdom_file.data_block.update(removed: true)

    expected = "Found 0 files to be archived.\n"

    expect {
      @cmd.execute
    }.to output(expected).to_stdout

    expect(wisdom_file.data_block.archive_id).to eq(nil)
  end
end
