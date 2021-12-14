require 'digest/md5'
require_relative '../lib/commands'

describe Metis::Assimilate do
  before(:each) do
    default_bucket('athena')
    @cmd = Metis::Assimilate.new
  end

  after(:each) do
    stubs.clear

    expect(stubs.contents(:athena)).to be_empty
  end

  it 'moves files and folders into the root folder path' do
    wisdom_path = stubs.create_data('stubs', 'wisdom.txt', WISDOM)
    helmet_path = stubs.create_data('stubs', 'blueprints/helmet.txt', HELMET)

    @cmd.execute('athena', 'files', '/', wisdom_path, ::File.dirname(helmet_path))

    expect(Metis::File.count).to eq(2)
    expect(Metis::Folder.count).to eq(1)

    helmet_file, wisdom_file = Metis::File.all.sort_by(&:file_name)
    blueprints_folder = Metis::Folder.first

    # the file records are there, nested appropriately, with real data
    expect(wisdom_file.file_name).to eq('wisdom.txt')
    expect(wisdom_file.folder).to be_nil
    expect(wisdom_file).to be_has_data
    expect(File.read(wisdom_file.data_block.location)).to eq(WISDOM)

    expect(helmet_file.file_name).to eq('helmet.txt')
    expect(helmet_file.folder).to eq(blueprints_folder)
    expect(helmet_file).to be_has_data
    expect(File.read(helmet_file.data_block.location)).to eq(HELMET)

    # folders are also created
    expect(blueprints_folder.folder_name).to eq('blueprints')
    expect(blueprints_folder).to be_root_folder

    # the original files are untouched
    expect(::File.exists?(wisdom_path)).to be_truthy
    expect(::File.exists?(helmet_path)).to be_truthy
    File.delete(wisdom_file.data_block.location)
    File.delete(helmet_file.data_block.location)
  end

  it 'moves files and folders into a folder path' do
    upload_folder = create_folder('athena', 'upload')
    stubs.create_folder('athena', 'files', 'upload')

    wisdom_path = stubs.create_data('stubs', 'wisdom.txt', WISDOM)
    helmet_path = stubs.create_data('stubs', 'blueprints/helmet.txt', HELMET)

    @cmd.execute('athena', 'files', '/upload', wisdom_path, ::File.dirname(helmet_path))

    expect(Metis::File.count).to eq(2)
    expect(Metis::Folder.count).to eq(2)

    upload_folder = Metis::Folder.where(folder_name: 'upload').first
    blueprints_folder = Metis::Folder.where(folder_name: 'blueprints').first
    wisdom_file = Metis::File.where(file_name: 'wisdom.txt').first
    helmet_file = Metis::File.where(file_name: 'helmet.txt').first

    expect(wisdom_file.folder).to eq(upload_folder)
    expect(helmet_file.folder).to eq(blueprints_folder)
    expect([wisdom_file, helmet_file]).to all(be_has_data)

    File.delete(wisdom_file.data_block.location)
    File.delete(helmet_file.data_block.location)
  end

  it 'skips over existing folders' do
    # there is a blueprints folder already
    blueprints_folder = create_folder('athena', 'blueprints')

    stubs.create_folder('athena', 'files', 'blueprints')

    # we try to assimilate an existing folder
    @cmd.execute('athena', 'files', '/', 'spec/data/stubs/blueprints')

    expect(Metis::File.count).to eq(0)
    expect(Metis::Folder.count).to eq(1)
  end

  it 'skips over existing files' do
    blueprints_folder = create_folder('athena', 'blueprints')
    stubs.create_folder('athena', 'files', 'blueprints')

    # file already exists in the blueprints folder
    helmet_file = create_file('athena', 'helmet.txt', HELMET, folder: blueprints_folder)
    stubs.create_file('athena', 'files', 'blueprints/helmet.txt', HELMET)

    helmet_path = stubs.create_data('stubs', 'helmet.txt', HELMET)

    # we try to assimilate the same file
    @cmd.execute('athena', 'files', '/blueprints', helmet_path)

    expect(Metis::File.count).to eq(1)
    expect(Metis::Folder.count).to eq(1)
  end

  it 'refuses to assimilate badly named files' do
    blueprints_folder = create_folder('athena', 'blueprints')
    stubs.create_folder('athena', 'files', 'blueprints')

    stubs.create_data('stubs', 'helmet::one.txt', HELMET)

    expect {
      @cmd.execute('athena', 'files', '/blueprints', 'spec/stubs/helmet::one.txt')
    }.to raise_error(ArgumentError)

    expect(Metis::Folder.count).to eq(1)
    expect(Metis::File.count).to eq(0)
  end
end
